""" Observational Data retrieval """
from __future__ import annotations

import os
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


def serialize_stations(
    stations: pd.DataFrame,
    path: os.PathLike[str],
    schism_station_flag: str = "1 0 0 0 0 0 0 0 0",
) -> None:
    """
    Serialize `stations` metadata to the provided `path`.

    The provided

    Example:
    --------

        >>> stations = pd.DataFrame({
        ...     'lon': [1., 11., 21.],
        ...     'lat': [1., 4., 1.],
        ...     'unique_id': ["a", "b", "c"],
        ...     'extra_col': ["AA", "BB", "CC"],
        ...     'mesh_index': [0, 1, 2],
        ...     'mesh_lon': [0., 10., 20.],
        ...     'mesh_lat': [0., 5., 0.],
        ...     'distance': [157249.38127194397, 157010.16264060183, 157249.38127194406],
        ... })
        >>> stations
           lon  lat unique_id extra_col  mesh_index  mesh_lon  mesh_lat       distance
        0    1    1         a        AA           0         0         0  157249.381272
        1   11    4         b        BB           1        10         5  157010.162641
        2   21    1         c        CC           2        20         0  157249.381272
        >>> serialize_stations(stations, "/tmp/station.in")

    Will create the following output:

        1 0 0 0 0 0 0 0 0	! https://schism-dev.github.io/schism/master/input-output/optional-inputs.html#stationin-bp-format
        3	! number of stations
        1 0.0000000000 0.0000000000 0 	!	0 1.0000000000 1.0000000000 157249.3812719440 AA a
        2 10.0000000000 5.0000000000 0 	!	1 11.0000000000 4.0000000000 157010.1626406018 BB b
        3 20.0000000000 0.0000000000 0 	!	2 21.0000000000 1.0000000000 157249.3812719441 CC c

    """
    # Sanity check
    mandatory_cols = {"mesh_index", "mesh_lon", "mesh_lat", "lon", "lat", "distance"}
    df_cols = set(stations.columns)
    if not df_cols.issuperset(mandatory_cols):
        msg = f"stations must have these columns too: {mandatory_cols.difference(df_cols)}"
        raise ValueError(msg)
    #
    basic_cols = ["mesh_lon", "mesh_lat", "z", "separator", "mesh_index", "lon", "lat", "distance"]
    extra_cols = df_cols.symmetric_difference(basic_cols).symmetric_difference(["z", "separator"])
    station_in_cols = basic_cols + list(sorted(extra_cols))
    station_in = stations.assign(
        z=0,
        separator="\t!\t",
    )
    station_in = station_in.set_index(station_in.index +1)
    station_in = station_in[station_in_cols]
    with open(f"{path}", "w") as fd:
        fd.write(f"{schism_station_flag.strip()}\t ! https://schism-dev.github.io/schism/master/input-output/optional-inputs.html#stationin-bp-format\n")
        fd.write(f"{len(station_in)}\t ! number of stations\n")
        station_in.to_csv(fd, header=None, sep=" ", float_format="%.10f")
