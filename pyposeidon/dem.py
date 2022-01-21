"""
Dem module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import importlib
import logging
import multiprocessing
import os
import pathlib
import warnings
import sys

import numpy as np
import pyresample
import xarray as xr

from pyposeidon.utils.fix import fix
from pyposeidon import tools

NCORES = max(1, multiprocessing.cpu_count() - 1)

logger = logging.getLogger(__name__)

LONGITUDE_NAMES = {"longitude", "lon", "x", "Lon", "LONGITUDE", "LON", "X"}
LATITUDE_NAMES = {"latitude", "lat", "y", "Lat", "LATITUDE", "LAT", "Y"}


class Dem:
    def __init__(self, dem_source=None, **kwargs):

        """
        Read bathymetric data from various sources.

        !!! danger ""
            Due to a limitation of the Library rendering the docstrings, all arguments are marked
            as `required`, nevertheless they are all `Optional` except geometry.

        Args:
            dem_source (str): Path or url to bathymetric data.
            lon_min (float) : Minimum longitude.
            lon_max (float) : Maximum longitude.
            lat_min (float) : Minimum latitude.
            lat_max (float) : Maximum latitude.
            geometry (Union[dict, str, GeoDataFrame]): A `GeoDataFrame` or the path to a shapefile or
                a dict defining the lat/lon window.
            coastlines (Union[str, GeoDataFrame]): A `GeoDataFrame` or the path to a shapefile which
                describes the coastlines. Defaults to `None`.
            adjust_dem (bool) :  Option to match dem with coastlines, Defaults to `True`.
            grid_x (list) : Array of longitude of mesh nodes.
            grid_y (list) : Array of latitude of mesh nodes.
            wet_only (bool) : Flag to mask dry values when interpolating on mesh. Defaults to `False`.
        """

        # integrate geometry attribute.
        geometry = kwargs.get("geometry", None)
        if isinstance(geometry, dict):
            kwargs.update(**geometry)

        if not dem_source:
            logger.error("dem_source is required\n")
            sys.exit(1)
        else:
            self.Dataset = dem_(source=dem_source, **kwargs)

        #        print('CHECK2  :','adjusted' in self.Dataset.data_vars.keys())
        #        print(self.Dataset.data_vars.keys())
        if kwargs.get("adjust_dem", True):
            coastline = kwargs.get("coastlines", None)
            if coastline is None:
                logger.warning("coastlines not present, aborting adjusting dem\n")
            elif "adjusted" in self.Dataset.data_vars.keys():
                logger.info("Dem already adjusted\n")
            elif "fval" in self.Dataset.data_vars.keys():
                logger.info("Dem already adjusted\n")
            else:
                self.adjust(coastline)

    def adjust(self, coastline, **kwargs):

        self.Dataset = fix(self.Dataset, coastline, **kwargs)
        try:
            check = np.isnan(self.Dataset.adjusted).sum() == 0
        except:
            check = np.isnan(self.Dataset.fval).sum() == 0

        if not check:
            logger.warning("Adjusting dem failed, keeping original values\n")
            self.Dataset = self.Dataset.drop_vars("adjusted")


def normalize_coord_names(dataset: xr.Dataset) -> xr.Dataset:
    """Return a dataset with coords containing "longitude" and "latitude" """
    coords = set(dataset.coords.keys())
    # longitude
    for lon_name in LONGITUDE_NAMES:
        if lon_name in coords:
            break
    else:
        raise ValueError(f"Couldn't normalize longitude: {coords}")
    # latitude
    for lat_name in LATITUDE_NAMES:
        if lat_name in coords:
            break
    else:
        raise ValueError(f"Couldn't normalize latitude: {coords}")
    dataset = dataset.rename({lon_name: "longitude", lat_name: "latitude"})
    return dataset


def normalize_elevation_name(dataset: xr.Dataset) -> xr.Dataset:
    data_vars = list(dataset.data_vars.keys())
    if (len(data_vars) > 1) & ("elevation" not in data_vars):
        warnings.warn(f"There are multiple data_vars. Assuming 'elevation' is the first one: {data_vars}")
    if "elevation" not in data_vars:
        dataset = dataset.rename({data_vars[0]: "elevation"})
    return dataset


def dem_(source=None, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, **kwargs) -> xr.Dataset:

    if isinstance(source, xr.Dataset):
        data = source
    else:
        dataset_kwargs = kwargs.pop("dem_xr_kwargs", {})
        data = tools.open_dataset(source, **dataset_kwargs)

    data = normalize_coord_names(data)
    data = normalize_elevation_name(data)

    dlon0 = round(data.longitude.data.min())
    dlon1 = round(data.longitude.data.max())

    # recenter the window
    if dlon1 - dlon0 == 360.0:

        lon0 = lon_min + 360.0 if lon_min < data.longitude.min() else lon_min
        lon1 = lon_max + 360.0 if lon_max < data.longitude.min() else lon_max

        lon0 = lon0 - 360.0 if lon0 > data.longitude.max() else lon0
        lon1 = lon1 - 360.0 if lon1 > data.longitude.max() else lon1

    else:

        lon0 = lon_min
        lon1 = lon_max

    if (lon_min < data.longitude.min()) or (lon_max > data.longitude.max()):
        logger.warning("Lon must be within {} and {}".format(data.longitude.min().values, data.longitude.max().values))
        logger.warning("compensating if global dataset available")

    if (lat_min < data.latitude.min()) or (lat_max > data.latitude.max()):
        logger.warning("Lat is within {} and {}".format(data.latitude.min().values, data.latitude.max().values))

    # get idx
    if lon_max - lon_min == dlon1 - dlon0:
        i0 = 0 if lon_min == dlon0 else int(data.longitude.shape[0] / 2) + 2  # compensate for below
        i1 = data.longitude.shape[0] if lon_max == dlon1 else -int(data.longitude.shape[0] / 2) - 2
    else:
        i0 = np.abs(data.longitude.data - lon0).argmin()
        i1 = np.abs(data.longitude.data - lon1).argmin()

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

    if i0 > i1:

        p1 = data.elevation.isel(longitude=slice(lon_0, data.longitude.size), latitude=slice(lat_0, lat_1))

        p1 = p1.assign_coords({"longitude": p1.longitude.values - 360.0})

        p2 = data.elevation.isel(longitude=slice(0, lon_1), latitude=slice(lat_0, lat_1))

        dem = xr.concat([p1, p2], dim="longitude")

    else:

        dem = data.elevation.isel(longitude=slice(lon_0, lon_1), latitude=slice(lat_0, lat_1))

    if np.abs(np.mean(dem.longitude) - np.mean([lon_min, lon_max])) > 170.0:
        c = np.sign(np.mean([lon_min, lon_max]))
        dem["longitude"] = dem["longitude"] + c * 360.0

    # ---------------------------------------------------------------------
    logger.info("dem done\n")
    # ---------------------------------------------------------------------

    dem_data = xr.merge([dem])

    grid_x = kwargs.get("grid_x", None)
    if grid_x is not None:
        dem_data = dem_on_mesh(dem_data, **kwargs)

    return dem_data


def dem_on_mesh(dataset, **kwargs):

    # ---------------------------------------------------------------------
    logger.info(".. interpolating on mesh ..\n")
    # ---------------------------------------------------------------------

    grid_x = kwargs.get("grid_x", None)
    grid_y = kwargs.get("grid_y", None)
    # resample on the given grid
    xx, yy = np.meshgrid(dataset.longitude, dataset.latitude)  # original grid

    # Translate for pyresample
    if xx.mean() < 0 and xx.min() < -180.0:
        xx = xx + 180.0
        gx = grid_x + 180.0
    elif xx.mean() > 0 and xx.max() > 180.0:
        xx = xx - 180.0
        gx = grid_x - 180.0
    else:
        gx = grid_x

    orig = pyresample.geometry.SwathDefinition(lons=xx, lats=yy)  # original points
    targ = pyresample.geometry.SwathDefinition(lons=gx, lats=grid_y)  # target grid

    var = "adjusted" if "adjusted" in dataset.data_vars.keys() else "elevation"

    # ---------------------------------------------------------------------
    logger.info(".. using {} values ..\n".format(var))
    # ---------------------------------------------------------------------

    wet = kwargs.get("wet_only", False)
    if wet:
        # mask positive bathymetry
        vals = np.ma.masked_array(dataset, dataset[var].values > 0)
    else:
        vals = dataset[var].values

    # with nearest using only the water values

    itopo = pyresample.kd_tree.resample_nearest(
        orig, dataset[var].values, targ, radius_of_influence=100000, fill_value=np.nan
    )  # ,nprocs=ncores)

    if len(grid_x.shape) > 1:
        idem = xr.Dataset(
            {
                "ival": (["k", "l"], itopo),
                "ilons": (["k", "l"], grid_x),
                "ilats": (["k", "l"], grid_y),
            }
        )  # ,
    #                           coords={'ilon': ('ilon', grid_x[0,:]),
    #                                   'ilat': ('ilat', grid_y[:,0])})

    elif len(grid_x.shape) == 1:
        idem = xr.Dataset(
            {
                "ival": (["k"], itopo),
                "ilons": (["k"], grid_x),
                "ilats": (["k"], grid_y),
            }
        )

    # ---------------------------------------------------------------------
    logger.info("dem on grid done\n")
    # ---------------------------------------------------------------------

    return xr.merge([dataset, idem])


def to_output(dataset, solver_name, **kwargs):
    solver_class = tools.get_solver(solver_name=solver_name)
    solver_class.to_dep(dataset, **kwargs)
