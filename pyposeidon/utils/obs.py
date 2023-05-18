"""
Observational Data retrieval

"""
from __future__ import annotations
from searvey import ioc
import pandas as pd
import geopandas as gp
from datetime import datetime


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

    return data
