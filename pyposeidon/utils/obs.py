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
    # merge
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
