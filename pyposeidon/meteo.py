"""
Meteo module. Pre-processing the weather forcing component.

"""

# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import datetime
import glob
import sys
import os
import pathlib
import time
from typing import Callable
import dask
import xarray as xr
import pandas as pd
import importlib
from pyposeidon.utils.get_value import get_value
from pyposeidon import tools
import logging

logger = logging.getLogger(__name__)


class Meteo:
    def __init__(self, meteo_source=None, **kwargs):

        """Read meteo data from variable sources.

        Parameters
        ----------
        meteo_source : list of str or url
            list of files, link
        engine : str
            Name of xarray backend to be used or 'url'

        Returns
        -------
        retrieved : xarray DataSet

        """

        # integrate geometry attribute.
        geometry = kwargs.get("geometry", None)
        if geometry:
            kwargs.update(**geometry)

        if isinstance(meteo_source, pathlib.Path):
            meteo_source = meteo_source.as_posix()

        meteo_func = dispatch_meteo_source(meteo_source)
        self.Dataset = meteo_func(meteo_source=meteo_source, **kwargs)

    def to_output(self, solver=None, **kwargs):

        model = importlib.import_module("pyposeidon.model")  # load pyposeidon model class

        s = getattr(model, solver)  # get solver class
        var_list = kwargs.pop("vars", ["msl", "u10", "v10"])

        m_index = get_value(self, kwargs, "m_index", 1)

        split_by = get_value(self, kwargs, "meteo_split_by", None)
        if split_by:
            times, datasets = zip(*self.Dataset.groupby("time.{}".format(split_by)))
            mpaths = ["sflux/sflux_air_{}.{:04d}.nc".format(m_index, t + 1) for t in np.arange(len(times))]
            for das, mpath in list(zip(datasets, mpaths)):
                s.to_force(das, vars=var_list, filename=mpath, **kwargs)
        else:
            s.to_force(self.Dataset, vars=var_list, **kwargs)


def passthrough(meteo_source, **kwargs):
    return meteo_source


def empty(meteo_source, **kwargs):
    return None


def cfgrib(
    meteo_source,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    start_date=None,
    end_date=None,
    time_frame=None,
    meteo_irange=[0, -1, 1],
    meteo_merge=None,
    meteo_combine_by="by_coords",
    **kwargs,
):

    backend_kwargs = kwargs.get("meteo_backend_kwargs", {"indexpath": ""})
    xr_kwargs = kwargs.get("meteo_xr_kwargs", {})  # {'concat_dim':'step'})
    minlon = lon_min
    maxlon = lon_max

    try:
        start_date = pd.to_datetime(start_date)
    except:
        pass

    if time_frame:
        try:
            end_date = start_date + pd.to_timedelta(time_frame)
        except:
            pass
    else:
        try:
            end_date = pd.to_datetime(end_date)
        except:
            pass

    ft1, ft2, dft = meteo_irange

    ts = pd.to_datetime(start_date)
    te = pd.to_datetime(end_date)

    # ---------------------------------------------------------------------
    logger.info("extracting meteo")
    # ---------------------------------------------------------------------

    data = xr.open_mfdataset(
        meteo_source,
        combine=meteo_combine_by,
        engine="cfgrib",
        backend_kwargs=backend_kwargs,
        **xr_kwargs,
    )

    data = data.squeeze(drop=True)
    #        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing> not efficient for output

    time_coord = data.msl.dims[0]

    if time_coord != "time":
        if "time" in data.coords.keys():
            data = data.rename({"time": "rtime"})
            data = data.rename({time_coord: "time"})
            data = data.assign_coords(time=data.valid_time)

    if meteo_merge == "last":
        logger.info("combining meteo by keeping the newest value")
        mask = data.time.to_pandas().duplicated("last").values
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            msl = data.msl[~mask]
            u10 = data.u10[~mask]
            v10 = data.v10[~mask]
        data = xr.merge([msl, u10, v10])

    elif meteo_merge == "first":
        logger.info("combining meteo by keeping the oldest value")
        mask = data.time.to_pandas().duplicated("last").values
        mask_ = np.array([mask[0]] + mask[:-1].tolist())
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            msl = data.msl[~mask_]
            u10 = data.u10[~mask_]
            v10 = data.v10[~mask_]
        data = xr.merge([msl, u10, v10])

    if not lon_min:
        lon_min = data.longitude.data.min()
    if not lon_max:
        lon_max = data.longitude.data.max()
    if not lat_min:
        lat_min = data.latitude.data.min()
    if not lat_max:
        lat_max = data.latitude.data.max()

    if lon_min < data.longitude.data.min():
        lon_min = lon_min + 360.0

    if lon_max < data.longitude.data.min():
        lon_max = lon_max + 360.0

    if lon_min > data.longitude.data.max():
        lon_min = lon_min - 360.0

    if lon_max > data.longitude.data.max():
        lon_max = lon_max - 360.0

    if not ts:
        ts = data.time.data[ft1]
    if not te:
        te = data.time.data[ft2]

    if ts < data.time.data.min():
        logger.warning(
            "coverage between {} and {} \n".format(
                pd.to_datetime(data.valid_time.values.min()).strftime("%Y.%m.%d %H:%M:%S"),
                pd.to_datetime(data.valid_time.values.max()).strftime("%Y.%m.%d %H:%M:%S"),
            )
        )
        logger.error("time frame does not match source range\n")
        sys.exit(1)

    if te > data.time.data.max():
        logger.warning(
            "coverage between {} and {} \n".format(
                pd.to_datetime(data.valid_time.values.min()).strftime("%Y.%m.%d %H:%M:%S"),
                pd.to_datetime(data.valid_time.values.max()).strftime("%Y.%m.%d %H:%M:%S"),
            )
        )
        logger.error("time frame does not match source range\n")
        sys.exit(1)

    tslice = slice(ts, te, dft)

    i0 = np.abs(data.longitude.data - lon_min).argmin()
    i1 = np.abs(data.longitude.data - lon_max).argmin()

    j0 = np.abs(data.latitude.data - lat_min).argmin()
    j1 = np.abs(data.latitude.data - lat_max).argmin()

    # expand the window a little bit
    lon_0 = max(0, i0 - 2)
    lon_1 = min(data.longitude.size, i1 + 2)

    lat_0 = max(0, j0 - 2)
    lat_1 = min(data.latitude.size, j1 + 2)

    # descenting lats
    if j0 > j1:
        j0, j1 = j1, j0
        lat_0 = max(0, j0 - 1)
        lat_1 = min(data.latitude.size, j1 + 3)

    if i0 >= i1:

        sh = (
            data[["msl", "u10", "v10"]]
            .isel(
                longitude=slice(lon_0, data.longitude.size),
                latitude=slice(lat_0, lat_1),
            )
            .sel(time=tslice)
        )
        sh = sh.assign_coords({"longitude": sh.longitude.values - 360.0})

        sh1 = (
            data[["msl", "u10", "v10"]].isel(longitude=slice(0, lon_1), latitude=slice(lat_0, lat_1)).sel(time=tslice)
        )

        tot = xr.concat([sh, sh1], dim="longitude")

    else:

        tot = (
            data[["msl", "u10", "v10"]]
            .isel(longitude=slice(lon_0, lon_1), latitude=slice(lat_0, lat_1))
            .sel(time=tslice)
        )

    # Adjust lon values
    try:
        if np.abs(np.mean([lon_min, lon_max]) - np.mean([minlon, maxlon])) > 300.0:
            c = np.sign(np.mean([minlon, maxlon]))
            tot["longitude"] = tot.longitude.values + c * 360.0
    except:
        pass

    # ---------------------------------------------------------------------
    logger.info("meteo done\n")
    # ---------------------------------------------------------------------

    return tot


def from_url(
    meteo_source,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    start_date=None,
    end_date=None,
    time_frame=None,
    **kwargs,
):

    if not start_date:
        start_date = pd.to_datetime(datetime.datetime.today())
    else:
        start_date = pd.to_datetime(start_date)

    if time_frame:
        end_date = start_date + pd.to_timedelta(time_frame)
    elif end_date:
        end_date = pd.to_datetime(end_date)

    ts = pd.to_datetime(start_date)
    te = pd.to_datetime(end_date)

    xr_kwargs = kwargs.get("meteo_xr_kwargs", {"engine": "pydap"})

    # ---------------------------------------------------------------------
    logger.info("extracting meteo from {}\n".format(meteo_source))
    # ---------------------------------------------------------------------

    data = xr.open_dataset(meteo_source, **xr_kwargs)

    try:
        data = data.rename({"lon": "longitude", "lat": "latitude"})
    except:
        pass

    lon0 = lon_min + 360.0 if lon_min < data.longitude.min() else lon_min
    lon1 = lon_max + 360.0 if lon_max < data.longitude.min() else lon_max

    lon0 = lon0 - 360.0 if lon0 > data.longitude.max() else lon0
    lon1 = lon1 - 360.0 if lon1 > data.longitude.max() else lon1

    # adjust te if value is None

    if not te:
        te = pd.to_datetime(data.time.max().values)

    if ts < data.time.min().values:
        ld = pd.to_datetime(data.time.min().values).strftime(format="%Y-%m-%d %H:%M:%S")
        hd = pd.to_datetime(data.time.max().values).strftime(format="%Y-%m-%d %H:%M:%S")
        logger.error("time frame not available\n")
        logger.warning("coverage between {} and {} \n".format(ld, hd))
        sys.exit(1)

    if te > data.time.max().values:
        logger.error("time frame not available\n")
        logger.warning("coverage between {} and {} \n".format(ld, hd))
        sys.exit(1)

    tslice = slice(ts, te)

    i0 = np.abs(data.longitude.data - lon0).argmin()
    i1 = np.abs(data.longitude.data - lon1).argmin()

    j0 = np.abs(data.latitude.data - lat_min).argmin()
    j1 = np.abs(data.latitude.data - lat_max).argmin()

    if i0 >= i1:

        sh = (
            data[["prmslmsl", "ugrd10m", "vgrd10m"]]
            .isel(longitude=slice(i0, data.longitude.size), latitude=slice(j0, j1 + 1))
            .sel(time=tslice)
        )
        sh = sh.assign_coords({"longitude": sh.longitude.values - 360.0})

        sh1 = (
            data[["prmslmsl", "ugrd10m", "vgrd10m"]]
            .isel(longitude=slice(0, i1 + 1), latitude=slice(j0, j1 + 1))
            .sel(time=tslice)
        )

        tot = xr.concat([sh, sh1], dim="longitude")

    else:

        tot = (
            data[["prmslmsl", "ugrd10m", "vgrd10m"]]
            .isel(longitude=slice(i0, i1 + 1), latitude=slice(j0, j1 + 1))
            .sel(time=tslice)
        )

    if np.abs(np.mean(tot.longitude) - np.mean([lon_min, lon_max])) > 300.0:
        c = np.sign(np.mean([lon_min, lon_max]))
        tot["longitude"] = tot["longitude"] + c * 360.0

    tot = tot.rename({"prmslmsl": "msl", "ugrd10m": "u10", "vgrd10m": "v10"})

    # ---------------------------------------------------------------------
    logger.info("meteo done\n")
    # ---------------------------------------------------------------------

    return tot


def netcdf(
    meteo_source,
    lon_min=None,
    lon_max=None,
    lat_min=None,
    lat_max=None,
    start_date=None,
    end_date=None,
    time_frame=None,
    meteo_irange=[0, -1, 1],
    meteo_combine_by="by_coords",
    **kwargs,
):

    # ---------------------------------------------------------------------
    logger.info("extracting meteo\n")
    # ---------------------------------------------------------------------

    xr_kwargs = kwargs.get("meteo_xr_kwargs", {})

    data = xr.open_mfdataset(meteo_source, combine=meteo_combine_by, **xr_kwargs)

    # rename var/coords
    time_coord = [x for x in data.coords if "time" in data[x].long_name.lower()]
    lon_coord = [x for x in data.coords if "longitude" in data[x].long_name.lower()]
    lat_coord = [x for x in data.coords if "latitude" in data[x].long_name.lower()]

    sn = {}
    for x in data.variables:
        try:
            sn.update({x: data[x].standard_name.lower()})
        except:
            sn.update({x: data[x].long_name.lower()})

    msl_ = [x for (x, v) in sn.items() if "pressure" in v]
    u10_ = [x for (x, v) in sn.items() if ("u wind" in v) | ("u-component" in v)]
    v10_ = [x for (x, v) in sn.items() if ("v wind" in v) | ("v-component" in v)]

    data = data.rename(
        {
            msl_[0]: "msl",
            u10_[0]: "u10",
            v10_[0]: "v10",
            lon_coord[0]: "longitude",
            lat_coord[0]: "latitude",
        }
    )

    data = data.sel(longitude=slice(lon_min, lon_max)).sel(latitude=slice(lat_min, lat_max))

    try:
        data = data.sel(time=slice(start_date, end_date))
    except:
        data = data.sel(time=slice(start_date, start_date + pd.to_timedelta(time_frame)))

    s, f, i = meteo_irange
    data = data.isel(time=slice(s, f, i))

    # ---------------------------------------------------------------------
    logger.info("meteo done\n")
    # ---------------------------------------------------------------------

    return data


def dispatch_meteo_source(obj) -> Callable[..., xr.Dataset]:
    """Choose the appropriate function to handle `obj`

    `obj` can be any of the following:

    - `str`
    - `Path`
    - `list[str]`
    - `list[Path]`
    - `None`
    - `URL`

    """
    original = obj
    obj = tools.cast_path_to_str(obj)
    if tools.is_iterable(obj) and not isinstance(obj, str) and not isinstance(obj, xr.Dataset):
        if len({type(o) for o in obj}) != 1:
            raise ValueError(f"Multiple types in iterable: {obj}")
        obj = tools.cast_path_to_str(obj[0])
    if obj is None:
        meteo_func = empty
    elif isinstance(obj, xr.Dataset):
        meteo_func = passthrough
    elif isinstance(obj, str):
        if obj.lower().endswith("grib"):
            meteo_func = cfgrib
        elif obj.lower().startswith("http"):
            meteo_func = from_url
        elif obj.lower().endswith("nc"):
            meteo_func = netcdf
        else:
            raise ValueError(f"Can't determine meteo_type from string: {obj}")
    else:
        raise ValueError(f"Can't determine meteo_type: {type(original)} - {original}")
    return meteo_func
