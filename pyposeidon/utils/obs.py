"""
Observational Data retrieval

"""
from __future__ import annotations
from searvey import ioc
import pandas as pd
import geopandas as gp
import numpy as np
import xarray as xr
from datetime import datetime
import itertools
from typing import Iterable
from typing import Iterator


def grouper(
    iterable: Iterable[_T],
    n: int,
    *,
    incomplete: str = "fill",
    fillvalue: Union[_U, None] = None,
) -> Iterator[Tuple[Union[_T, _U], ...]]:
    """Collect data into non-overlapping fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == "fill":
        return itertools.zip_longest(*args, fillvalue=fillvalue)
    if incomplete == "strict":
        return zip(*args, strict=True)  # type: ignore[call-overload]
    if incomplete == "ignore":
        return zip(*args)
    else:
        raise ValueError("Expected fill, strict, or ignore")


def merge_datasets(datasets: List[xr.Dataset], size: int = 5) -> List[xr.Dataset]:
    datasets = [xr.merge(g for g in group if g) for group in grouper(datasets, size)]
    return datasets


def get_bogus(temp, idfs, ntime):
    return xr.Dataset(
        {"nodata": (["ioc_code", "time"], temp)},
        coords={
            "ioc_code": (["ioc_code"], [idfs]),
            "time": ntime,
        },
    )


def get_obs_data(stations: str | gp.GeoDataFrame, start_time=None, end_time=None, period=None, **kwargs):
    if isinstance(stations, str):
        df = gp.read_file(stations)
    elif isinstance(stations, gp.GeoDataFrame):
        df = stations
    else:
        raise ValueError("No valid DataFrame")

    if not period:
        dt = pd.to_datetime(end_time) - pd.to_datetime(start_time)
        period = dt.total_seconds() / 3600 / 24

    data = ioc.get_ioc_data(
        ioc_metadata=df,
        endtime=datetime(end_time.year, end_time.month, end_time.day, end_time.hour),
        period=period,
    )

    if "id" not in stations.columns:
        stations["id"] = [f"IOC-{x}" for x in stations.ioc_code]

    s2 = stations.loc[stations.ioc_code.isin(data.ioc_code.values)].id.values

    # create bogus data when no data
    sn = stations.loc[~stations.ioc_code.isin(s2)]
    sa = sn[["ioc_code", "lat", "lon", "country", "location"]].to_xarray()
    ntime = data.isel(ioc_code=0).time.data
    temp = np.array([np.NaN] * ntime.shape[0])[np.newaxis, :]  # nans for the time frame
    xdfs = []
    for idfs in sn.ioc_code:
        xdfs.append(get_bogus(temp, idfs, ntime))

    # in order to keep memory consumption low, let's group the datasets
    # and merge them in batches
    while len(xdfs) > 5:
        xdfs = merge_datasets(xdfs)
    # Do the final merging
    sns = xr.merge(xdfs)

    # set correct dims for metadata
    sa = sa.set_coords("ioc_code").swap_dims({"index": "ioc_code"}).reset_coords("index").drop_vars("index")

    # merge nodata with metadata
    bg = xr.merge([sns, sa])

    # merge bogus with searvey data
    data = xr.merge([data, bg])

    sts = stations.sort_values(by=["ioc_code"])  # make sure the indexing works

    data = data.assign_coords({"id": ("ioc_code", sts.id)}).swap_dims({"ioc_code": "id"}).reset_coords("ioc_code")

    return data
