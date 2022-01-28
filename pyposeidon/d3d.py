"""
Main d3d module of pyposeidon. It controls the creation, output & execution of a complete simulation based on DELFT3D

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import os
import datetime
import numpy as np
import xml.dom.minidom as md
from shutil import copy2
import subprocess
import sys
import json
from collections import OrderedDict
import pandas as pd
import glob
from shutil import copyfile
import xarray as xr
import geopandas as gp

# local modules
import pyposeidon
import pyposeidon.mesh as pmesh
import pyposeidon.meteo as pmeteo
import pyposeidon.dem as pdem
from pyposeidon.paths import DATA_PATH
from pyposeidon.utils.get_value import get_value
from pyposeidon.utils.converter import myconverter
from pyposeidon.utils import data
import logging

from .bnd import Box

logger = logging.getLogger(__name__)

import multiprocessing

NCORES = max(1, multiprocessing.cpu_count() - 1)

# strings to be used
le = ["A", "B"]

nm = ["Z", "A"]


D3D_NAME = "d3d"


class d3d:
    def __init__(self, **kwargs):

        """
        Create a D3D solver

        !!! danger ""
            Due to a limitation of the Library rendering the docstrings, all arguments are marked
            as `required`, nevertheless they are all `Optional`.

        Args:
            geometry Union[dict, str, GeoDataFrame]: A `GeoDataFrame` or the path to a shapefile or
                a dict defining the lat/lon window.
            start_date str: The date from which the analysis should start. It should be a string parseable
                by `pd.to_datetime()`.
            end_date str: The date at which the analysis should end. It should be a string parseable by
                `pd.to_datetime()`.
            time_frame str: The duration of the analysis. It should be a string parseable by
                `pd.to_datetime()`.
            date str: Reference date of the run.
            meteo_source str: Path or url to meteo data.
            dem_source str: Path or url to bathymetric data.
            argfile str: Path to `_hydro.xml` file.
            update str: Control the update of the model e.g `['dem']`-> updates only bathymetry.
                Defaults to `["all"]`.
            rpath str: Path for output of the model. Defaults to `./d3d/`.
            tide str: Flag indicating whether to load "tide". Defaults to `False`.
            atm bool: The solver's atm. Defaults to `True`.
            tag str: The model's "tag". Defaults to `"d3d"`.
            resolution float: size of the regular grid. Defaults to `0.1`.
            ofilename str: Path to station file. Defaults to `None`.
            epath str: The path to the schism executable. If the `D3D` env variable has been
                set, then it overrides the value passed as the parameter.
            config_file str: Path to mdf file. Defaults to `None`.
            config dict: Parameter options passed to mdf file.
            output bool: Flag for saving to file. Defaults to `False`.
            update list[str]: Control the update of the model e.g `['dem']`-> updates only bathymetry.
                Defaults to `["all"]`.
        """

        self.geometry = kwargs.get("geometry", None)

        if self.geometry:

            if isinstance(self.geometry, dict):
                self.lon_min = self.geometry["lon_min"]
                self.lon_max = self.geometry["lon_max"]
                self.lat_min = self.geometry["lat_min"]
                self.lat_max = self.geometry["lat_max"]
            elif isinstance(self.geometry, str):

                try:
                    geo = gp.GeoDataFrame.from_file(self.geometry)
                except:
                    logger.error("geometry argument not a valid geopandas file")
                    sys.exit(1)

                (
                    self.lon_min,
                    self.lat_min,
                    self.lon_max,
                    self.lat_max,
                ) = geo.total_bounds

        start_date = kwargs.get("start_date", None)
        self.start_date = pd.to_datetime(start_date)

        if "time_frame" in kwargs:
            time_frame = kwargs.get("time_frame", None)
            self.end_date = self.start_date + pd.to_timedelta(time_frame)
            self.time_frame = time_frame
        elif "end_date" in kwargs:
            end_date = kwargs.get("end_date", None)
            self.end_date = pd.to_datetime(end_date)
            self.time_frame = self.end_date - self.start_date

        if not hasattr(self, "date"):
            self.date = self.start_date

        if not hasattr(self, "end_date"):
            # ---------------------------------------------------------------------
            logger.warning("model not set properly, No end_date\n")
            # ---------------------------------------------------------------------

        self.tag = kwargs.get("tag", "d3d")
        self.resolution = kwargs.get("resolution", 0.1)
        self.irange = kwargs.get("irange", [0, -1, 1])
        self.tide = kwargs.get("tide", False)
        self.atm = kwargs.get("atm", True)
        self.ofilename = kwargs.get("ofilename", None)

        self.solver_name = D3D_NAME

        try:
            self.epath = os.environ["D3D"]
        except:
            self.epath = kwargs.get("epath", None)

        for attr, value in kwargs.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

    # ============================================================================================
    # CONFIG
    # ============================================================================================
    def config(self, **kwargs):

        mdf_file = kwargs.get("config_file", None)
        dic = get_value(self, kwargs, "parameters", None)

        if mdf_file:
            self.mdf = pd.read_csv(mdf_file, sep="=")
        else:
            self.mdf = pd.read_csv(DATA_PATH + "default.mdf", sep="=")

        self.mdf = self.mdf.set_index(self.mdf.columns[0])  # set index

        mdfidx = self.mdf.index.str.strip()  # store the stripped names

        # define mesh file
        self.mdf.loc[self.mdf.index.str.contains("Filcco")] = "#{}#".format(self.tag + ".grd")

        # define enc file
        self.mdf.loc[self.mdf.index.str.contains("Filgrd")] = "#{}#".format(self.tag + ".enc")

        # define dep file
        self.mdf.loc[self.mdf.index.str.contains("Fildep")] = "#{}#".format(self.tag + ".dep")

        # define obs file
        if self.ofilename:
            self.mdf.loc[self.mdf.index.str.contains("Filsta")] = "#{}#".format(self.tag + ".obs")
        else:
            self.mdf.loc[self.mdf.index.str.contains("Filsta")] = "##"

        # adjust ni,nj
        nj, ni = self.nj, self.ni
        self.mdf.loc[self.mdf.index.str.contains("MNKmax")] = "{} {} {}".format(ni + 1, nj + 1, 1)  # add one like ddb

        # adjust iteration date
        self.mdf.loc[self.mdf.index.str.contains("Itdate")] = "#{}#".format(self.date.strftime(format="%Y-%m-%d"))

        # set time unit
        self.mdf.loc[self.mdf.index.str.contains("Tunit")] = "#M#"

        # adjust iteration start
        Tstart = self.start_date.hour * 60
        self.mdf.loc[self.mdf.index.str.contains("Tstart")] = Tstart

        # adjust iteration stop
        Tstop = Tstart + int(pd.to_timedelta(self.time_frame).total_seconds() / 60)
        self.mdf.loc[self.mdf.index.str.contains("Tstop")] = Tstop

        # adjust time for output
        mstep = get_value(self, kwargs, "map_step", 60)
        hstep = get_value(self, kwargs, "his_step", 0)
        pstep = get_value(self, kwargs, "pp_step", 0)
        rstep = get_value(self, kwargs, "restart_step", 0)

        if rstep == -1:  # save a restart file at the end
            rstep = Tstop

        self.mdf.loc[self.mdf.index.str.contains("Flmap")] = "{:d} {:d} {:d}".format(Tstart, mstep, Tstop)
        self.mdf.loc[self.mdf.index.str.contains("Flhis")] = "{:d} {:d} {:d}".format(Tstart, hstep, Tstop)
        self.mdf.loc[self.mdf.index.str.contains("Flpp")] = "{:d} {:d} {:d}".format(Tstart, pstep, Tstop)
        self.mdf.loc[self.mdf.index.str.contains("Flrst")] = rstep

        # time interval to smooth the hydrodynamic boundary conditions
        self.mdf.loc[self.mdf.index.str.contains("Tlfsmo")] = 0.0

        if not self.atm:
            self.mdf.loc["Sub1"] = " "

        # set tide only run
        if self.tide:
            self.mdf.loc[self.mdf.index.str.contains("Filbnd")] = "#{}#".format(self.tag + ".bnd")
            self.mdf.loc[self.mdf.index.str.contains("Filana")] = "#{}#".format(self.tag + ".bca")
        #           if 'Tidfor' not in order: order.append('Tidfor')
        #           inp['Tidfor']=[['M2','S2','N2','K2'], \
        #                       ['K1','O1','P1','Q1'], \
        #                         ['-----------']]

        # specify ini file
        # if 'Filic' not in order: order.append('Filic')
        # inp['Filic']=basename+'.ini'

        # netCDF output
        if not "FlNcdf" in mdfidx:
            self.mdf.reindex(self.mdf.index.values.tolist() + ["FlNcdf "])

        self.mdf.loc["FlNcdf "] = "#map his#"

        other = kwargs.get("config", None)
        if other:
            # Check for any other mdf variable in input
            for key, val in other.items():
                if key in mdfidx:
                    self.mdf.loc[self.mdf.index.str.contains(key)] = val
                else:
                    self.mdf.loc[key] = val

        output = kwargs.get("output", False)

        if output:
            # save mdf
            path = get_value(self, kwargs, "rpath", "./d3d/")
            self.mdf.to_csv(path + self.tag + ".mdf", sep="=")

    # ============================================================================================
    # METEO
    # ============================================================================================

    def force(self, **kwargs):

        meteo_source = get_value(self, kwargs, "meteo_source", None)

        kwargs.update({"meteo_source": meteo_source})

        flag = get_value(self, kwargs, "update", [])
        # check if files exist

        z = {**self.__dict__, **kwargs}  # merge self and possible kwargs

        if flag:
            if ("meteo" in flag) | ("all" in flag):
                self.meteo = pmeteo.Meteo(**z)
            else:
                logger.info("skipping meteo files ..\n")
        else:
            self.meteo = pmeteo.Meteo(**z)

    @staticmethod
    def from_force(filename=None, name=None):

        df = pd.read_csv(filename, header=0, names=["data"], index_col=None, low_memory=False)

        tlines = df[df.data.str.contains("TIME")].index  # rows which start with TIME

        # get attrs
        d1 = df.loc[0 : tlines[0] - 1, "data"].str.split("=", 2, expand=True)
        d1.columns = ["key", "value"]  # assign column names
        d1.key = d1.key.str.strip()  # cleanup spaces
        d1.value = d1.value.str.strip()
        attrs = dict(zip(d1.key, d1.value))  # create dict
        for key in ["n_cols", "n_rows", "n_quantity"]:  # str -> int
            attrs[key] = int(attrs[key])

        for key in ["x_llcenter", "dx", "y_llcenter", "dy", "NODATA_value"]:
            attrs[key] = float(attrs[key])

        # get time reference
        d2 = df.loc[tlines, "data"].str.split("=", 2, expand=True)
        d2 = d2.drop(d2.columns[0], axis=1)
        d2.columns = ["data"]
        d2 = d2.loc[:, "data"].str.split(" ", 4, expand=True)
        d2 = d2.drop(d2.columns[[0, 2, 3]], axis=1)
        d2.columns = ["hours", "time0"]
        d2.hours = d2.hours.apply(pd.to_numeric)
        d2.time0 = pd.to_datetime(d2.time0.values)
        d2 = d2.reset_index(drop=True)
        # create timestamps
        time = []
        for i in range(d2.shape[0]):
            time.append(d2.time0[0] + pd.DateOffset(hours=int(d2.loc[i, "hours"])))
        d2["time"] = time

        # get the float numbers
        d3 = df.drop(np.arange(0, tlines[0]))
        d3 = d3.drop(tlines)

        #    data = []
        #    for i in range(d3.values.shape[0]):
        #        row = d3.values[i][0].split(' ')
        #        row = [np.float(x) for x in row]
        #        data.append(row)
        #    data = np.array(data) # make array

        data = d3[d3.columns[0]].str.split(" ", attrs["n_cols"], expand=True).to_numpy().astype(float)

        data = data.reshape(d2.shape[0], attrs["n_rows"], attrs["n_cols"])  # reshape

        # define lat/lon
        lon = [attrs["x_llcenter"] + attrs["dx"] * i for i in np.arange(attrs["n_cols"])]
        lat = [attrs["y_llcenter"] + attrs["dy"] * i for i in np.arange(attrs["n_rows"])]

        # create an xarray
        da = xr.DataArray(
            data,
            dims=["time", "latitude", "longitude"],
            coords={"time": d2.time, "latitude": lat, "longitude": lon},
            name=name,
        )

        da.attrs = attrs

        return da

    @staticmethod
    def to_force(ar, **kwargs):

        logger.info("writing meteo files ..\n")

        path = kwargs.get("rpath", "./d3d/")

        [p, u, v] = kwargs.get("vars", "[None,None,None]")

        curvi = kwargs.get("curvi", False)

        flip = np.diff(ar.latitude.values)[0]

        dlat = np.abs(flip)
        dlon = np.diff(ar.longitude.values)[0]
        lat0 = ar.latitude.data.min()
        lon0 = ar.longitude.data.min()

        nodata = -9999.000

        pp = ar[p].fillna(nodata).values
        uu = ar[u].fillna(nodata).values
        vv = ar[v].fillna(nodata).values

        if not os.path.exists(path):
            os.makedirs(path)

            # open files
        pfid = open(path + "p.amp", "w")
        ufid = open(path + "u.amu", "w")
        vfid = open(path + "v.amv", "w")

        fi = [pfid, ufid, vfid]
        wi = [ufid, vfid]

        # write file headers
        for f in fi:
            f.write("FileVersion      = 1.03\n")
        if curvi:
            for f in fi:
                f.write("Filetype         = meteo_on_curvilinear_grid\n")
                f.write("grid_file        = wind.grd\n")
                f.write("first_data_value = grid_ulcorner\n")
                f.write("data_row         = grid_row\n")
        else:
            for f in fi:
                f.write("Filetype         = meteo_on_equidistant_grid\n")
                f.write("n_cols           = {}\n".format(ar[u].shape[2]))
                f.write("n_rows           = {}\n".format(ar[u].shape[1]))
                f.write("grid_unit        = degree\n")
                # code currently assumes lon and lat are increasing
                f.write("x_llcenter       = {:g}\n".format(lon0))
                f.write("dx               = {:g}\n".format(dlon))
                f.write("y_llcenter       = {:g}\n".format(lat0))
                f.write("dy               = {:g}\n".format(dlat))

        for f in fi:
            f.write("NODATA_value     = {:.3f}\n".format(nodata))
            f.write("n_quantity       = 1\n")

        ufid.write("quantity1        = x_wind\n")
        vfid.write("quantity1        = y_wind\n")
        pfid.write("quantity1        = air_pressure\n")

        for f in wi:
            f.write("unit1            = m s-1\n")

        pfid.write("unit1            = Pa\n")

        time0 = pd.to_datetime("2000-01-01 00:00:00")

        # write time blocks
        indx = ar.time.values - time0.to_datetime64()
        indx = indx.astype("timedelta64[m]") / 60

        for it in range(indx.size):  # nt + 0 hour
            for f in fi:
                f.write("TIME = {} hours since 2000-01-01 00:00:00 +00:00\n".format(indx[it].astype(int)))

            if flip < 0:
                np.savetxt(pfid, np.flipud(pp[it, :, :]), fmt="%.8f")
                np.savetxt(ufid, np.flipud(uu[it, :, :]), fmt="%.8f")
                np.savetxt(vfid, np.flipud(vv[it, :, :]), fmt="%.8f")
            else:
                np.savetxt(pfid, pp[it, :, :], fmt="%.8f")
                np.savetxt(ufid, uu[it, :, :], fmt="%.8f")
                np.savetxt(vfid, vv[it, :, :], fmt="%.8f")

        # close files
        for f in fi:
            f.close()

    # ============================================================================================
    # DEM
    # ============================================================================================
    @staticmethod
    def from_dep(filename, **kwargs):

        rdem = np.loadtxt(filename)

        dr = xr.DataArray(rdem[:-1, :-1], name="ival", dims=["k", "l"])

        return dr

    def bath(self, **kwargs):

        kwargs["grid_x"] = self.mesh.Dataset.lons.values
        kwargs["grid_y"] = self.mesh.Dataset.lats.values

        dpath = get_value(self, kwargs, "dem_source", None)

        kwargs.update({"dem_source": dpath})

        flag = get_value(self, kwargs, "update", [])
        # check if files exist
        if flag:
            if ("dem" in flag) | ("all" in flag):
                kwargs.update(
                    {
                        "lon_min": self.lon_min,
                        "lat_min": self.lat_min,
                        "lon_max": self.lon_max,
                        "lat_max": self.lat_max,
                    }
                )
                self.dem = pdem.Dem(**kwargs)
            else:
                logger.info("reading local dem file ..\n")
                dem_source = z["rpath"] + self.tag + ".dep"
                rdem = from_dep(dem_source)

        else:
            kwargs.update(
                {
                    "lon_min": self.lon_min,
                    "lat_min": self.lat_min,
                    "lon_max": self.lon_max,
                    "lat_max": self.lat_max,
                }
            )
            self.dem = pdem.Dem(**kwargs)

    @staticmethod
    def to_dep(dr, dry_mask=True, **kwargs):
        # save dem
        logger.info("writing dem file ..\n")
        path = kwargs.get("rpath", "./d3d/")

        flag = kwargs.get("update", None)
        tag = kwargs.get("tag", "d3d")

        try:
            try:
                bat = -dr.fval.values.astype(float)  # reverse for the hydro run/use the adjusted values
            #     mask = bat==999999
            except AttributeError:
                bat = -dr.ival.values.astype(float)  # reverse for the hydro run/revert to interpolated values

            nj, ni = bat.shape

            if dry_mask:

                mask = ~np.isnan(bat)  # mask out potential nan points
                mask[mask] = np.less(bat[mask], 0)  # get mask for dry points

                bat[mask] = np.nan  # mask dry points

            # append the line/column of nodata
            nodata = np.empty(ni)
            nodata.fill(np.nan)
            bat1 = np.vstack((bat, nodata))
            nodata = np.empty((nj + 1, 1))
            nodata.fill(np.nan)
            bat2 = np.hstack((bat1, nodata))

            bat2[np.isnan(bat2)] = -999.0

        except AttributeError:
            logger.warning("problem with dem Dataset ..")

        # Write bathymetry file
        if flag:
            if ("all" in flag) or ("dem" in flag):
                np.savetxt(path + tag + ".dep", bat2)
            else:
                logger.info("keeping dem file ..\n")
        else:
            np.savetxt(path + tag + ".dep", bat2)

    # ============================================================================================
    # BOUNDARY CONDITIONS TODO
    # ============================================================================================

    def bc(self, **kwargs):
        # define boundaries
        z = self.__dict__.copy()

        z["lons"] = self.mesh.Dataset.lons[0, :]
        z["lats"] = self.mesh.Dataset.lats[:, 0]

        try:
            ba = -self.dem.Dataset.ival.astype(float)
            # ba[ba<0]=np.nan
            z["dem"] = ba
            z["cn"] = 10

            z.update(kwargs)

            self.bound = Box(**z)

        except:
            logger.info("boundary files not set..\n")

    def to_bnd(self):
        # save bnd
        with open(path + self.tag + ".bnd", "w") as f:

            dd = OrderedDict(
                [
                    ("North", self.bound.North),
                    ("South", self.bound.South),
                    ("West", self.bound.West),
                    ("East", self.bound.East),
                ]
            )

            #    for key,val in self.bound.__dict__.items():
            for i, (key, val) in enumerate(dd.items()):  # to match deltares

                idx = 1
                for k1, k2 in val:
                    bname = key + str(idx)
                    f.write(
                        "{0:<10s}{1:>12s}{2:>2s}{3:>6d}{4:>6d}{5:>6d}{6:>6d}   0.0000000e+00 {7:<s}{8:<g}A {9:<s}{10:<g}B\n".format(
                            bname,
                            nm[0],
                            nm[1],
                            k1[0] + 1,
                            k1[1] + 1,
                            k2[0] + 1,
                            k2[1] + 1,
                            key,
                            idx,
                            key,
                            idx,
                        )
                    )  # fortran index ??
                    idx += 1

    def to_bca(self):
        # save bca
        with open(path + self.tag + ".bca", "w") as f:

            dd = OrderedDict(
                [
                    ("North", self.tide.North),
                    ("South", self.tide.South),
                    ("West", self.tide.West),
                    ("East", self.tide.East),
                ]
            )

            #     for key,val in self.tide.__dict__.items():
            for i, (key, val) in enumerate(dd.items()):  # to match deltares

                idx = 1
                if val:
                    l = np.arange(val.ampl.shape[0]) + idx
                    nl = [x for pair in zip(l, l) for x in pair]
                    sl = val.ampl.shape[0] * le
                    for t1, t2, amp, phase in zip(np.transpose(nl), np.transpose(sl), val.ampl, val.phase):
                        f.write("{}{}{}\n".format(key, t1, t2))
                        for a, b, c in zip(val.constituents, amp.flatten(), phase.flatten()):
                            f.write("{0:<3s}        {1:<.7e}   {2:<.7e}\n".format(a, b, c))

    def tidebc(self, **kwargs):

        self.tide = tide()
        for key, val in self.bound.__dict__.items():

            # compute tide constituents
            tval = []
            if len(val) > 0.0:
                blons = []
                blats = []
                for l1, l2 in val:
                    blons.append(self.mesh.Dataset.lons[l1[1] - 1, l1[0] - 1])
                    blats.append(self.mesh.Dataset.lats[l1[1] - 1, l1[0] - 1])
                    blons.append(self.mesh.Dataset.lons[l2[1] - 1, l2[0] - 1])
                    blats.append(self.mesh.Dataset.lats[l2[1] - 1, l2[0] - 1])

                blons = np.array(blons)  # .ravel().reshape(-1,2)[:,0]
                blats = np.array(blats)  # .ravel().reshape(-1,2)[:,1]
                #                  print(bound,blons,blats)

                tval = tide(tmodel=self.tmodel, tpath=self.tpath, blons=blons, blats=blats)

            setattr(self.tide, key, tval)

    @staticmethod
    def to_obs(self, **kwargs):
        # save obs

        ofilename = get_value(self, kwargs, "ofilename", None)
        flag = get_value(self, kwargs, "update", [])

        if ofilename:

            obs_points = pd.read_csv(
                ofilename,
                delimiter="\t",
                header=None,
                names=["index", "Name", "lat", "lon"],
            )
            obs_points = obs_points.set_index("index", drop=True).reset_index(drop=True)  # reset index if any

            obs_points = obs_points[
                (
                    obs_points.lon.between(
                        self.mesh.Dataset.lons.values.min(),
                        self.mesh.Dataset.lons.values.max(),
                    )
                )
                & (
                    obs_points.lat.between(
                        self.mesh.Dataset.lats.values.min(),
                        self.mesh.Dataset.lats.values.max(),
                    )
                )
            ]

            obs_points.reset_index(inplace=True, drop=True)

            try:
                bat = -self.dem.Dataset.fval.values.astype(float)  # reverse for the hydro run/use the adjusted values
            #     mask = bat==999999
            except AttributeError:
                bat = -self.dem.Dataset.ival.values.astype(
                    float
                )  # reverse for the hydro run/revert to interpolated values

            b = np.ma.masked_array(bat, np.isnan(bat))  # mask land

            i_indx, j_indx = self.vpoints(self.mesh.Dataset, obs_points, b, **kwargs)

            obs_points["i"] = i_indx
            obs_points["j"] = j_indx

            # drop NaN points
            obs = obs_points.dropna().copy()

            obs = obs.reset_index(drop=True)  # reset index

            obs["i"] = obs["i"].values.astype(int)
            obs["j"] = obs["j"].values.astype(int)
            obs["new_lat"] = self.mesh.Dataset.y[obs.i.values].values  # Valid point
            obs["new_lon"] = self.mesh.Dataset.x[obs.j.values].values

            self.obs = obs  # store it

            obs.Name = obs.Name.str.strip().apply(lambda name: name.replace(" ", ""))  # Remove spaces to write to file
            sort = sorted(obs.Name.values, key=len)  # sort the names to get the biggest word
            try:
                wsize = len(sort[-1])  # size of bigget word in order to align below
            except:
                pass

        if flag:

            if ("all" in flag) | ("model" in flag):

                # Add one in the indices due to python/fortran convention
                try:
                    with open(self.rpath + "{}.obs".format(self.tag), "w") as f:
                        for l in range(obs.shape[0]):
                            f.write(
                                "{0:<{3}}{1:>{3}}{2:>{3}}\n".format(
                                    obs.Name[l][:20], obs.j[l] + 1, obs.i[l] + 1, wsize
                                )
                            )
                except:  # TODO
                    pass

        else:
            try:
                # Add one in the indices due to python/fortran convention
                with open(self.rpath + "{}.obs".format(self.tag), "w") as f:
                    for l in range(obs.shape[0]):
                        f.write(
                            "{0:<{3}}{1:>{3}}{2:>{3}}\n".format(obs.Name[l][:20], obs.j[l] + 1, obs.i[l] + 1, wsize)
                        )
            except:
                pass

    # ============================================================================================
    # EXECUTION
    # ============================================================================================
    def create(self, **kwargs):

        if not kwargs:
            kwargs = self.__dict__.copy()

        # Grid
        self.mesh = pmesh.set(type="r2d", **kwargs)

        # set lat/lon from file
        if hasattr(self, "mesh_file"):
            kwargs.update({"lon_min": self.mesh.Dataset.x.values.min()})
            kwargs.update({"lon_max": self.mesh.Dataset.x.values.max()})
            kwargs.update({"lat_min": self.mesh.Dataset.y.values.min()})
            kwargs.update({"lat_max": self.mesh.Dataset.y.values.max()})

        nj, ni = self.mesh.Dataset.lons.shape
        self.nj, self.ni = nj, ni

        kwargs.update({"ni": ni, "nj": nj})

        # get bathymetry
        self.bath(**kwargs)

        # get boundaries
        self.bc()

        # get meteo
        if self.atm:
            self.force(**kwargs)

        # get tide
        if self.tide:
            self.tidebc()

        self.config(**kwargs)

    def run(self, **kwargs):

        calc_dir = get_value(self, kwargs, "rpath", "./d3d/")

        try:
            bin_path = os.environ["D3D"]
        except:
            bin_path = get_value(self, kwargs, "epath", None)

        try:
            lib_path = os.environ["LD3D"]
        except:
            lib_path = get_value(self, kwargs, "lpath", None)

        if bin_path is None:
            # ------------------------------------------------------------------------------
            logger.warning("D3D executable path (epath) not given -> using default \n")
            # ------------------------------------------------------------------------------
            bin_path = os.pathsep + cpath
            lib_path = bin_path

        ncores = get_value(self, kwargs, "ncores", NCORES)

        argfile = get_value(self, kwargs, "argfile", self.tag + "_hydro.xml")

        # ---------------------------------------------------------------------
        logger.info("executing model\n")
        # ---------------------------------------------------------------------

        # note that cwd is the folder where the executable is
        ex = subprocess.Popen(
            args=["./run_flow2d3d.sh {} {} {}".format(argfile, ncores, bin_path, lib_path)],
            cwd=calc_dir,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )  # , bufsize=1)

        with open(calc_dir + self.tag + "_run.log", "w") as f:  # save output

            for line in iter(ex.stdout.readline, b""):
                f.write(line.decode(sys.stdout.encoding))
            #                logger.info(line.decode(sys.stdout.encoding))

            for line in iter(ex.stderr.readline, b""):
                logger.info(line.decode(sys.stdout.encoding))
                tempfiles = glob.glob(calc_dir + "/tri-diag." + self.tag + "-*")
                try:
                    biggest = max(tempfiles, key=(lambda tf: os.path.getsize(tf)))
                    with open(biggest, "r") as f1:
                        for line in f1:
                            f.write(line.decode(sys.stdout.encoding))
                except:
                    pass

        # cleanup
        tempfiles = glob.glob(calc_dir + "/tri-diag." + self.tag + "-*")
        biggest = max(tempfiles, key=(lambda tf: os.path.getsize(tf)))
        with open(calc_dir + self.tag + "_run.log", "a") as f:  # save diagnosis
            with open(biggest, "r") as f1:
                for line in f1:
                    f.write(line)

        tempfiles = glob.glob(calc_dir + "/tri-diag." + self.tag + "-*") + glob.glob(calc_dir + "/TMP_*")

        for filename in tempfiles:
            try:
                os.remove(filename)
            except OSError:
                pass

        ex.stdout.close()
        ex.stderr.close()

        # ---------------------------------------------------------------------
        logger.info("FINISHED\n")
        # ---------------------------------------------------------------------

    def save(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./d3d/")

        lista = [key for key, value in self.__dict__.items() if key not in ["meteo", "dem", "mesh"]]
        dic = {k: self.__dict__.get(k, None) for k in lista}

        mesh = self.__dict__.get("mesh", None)
        if isinstance(mesh, str):
            dic.update({"mesh": mesh})
        else:
            dic.update({"mesh": mesh.__class__.__name__})

        dem = self.__dict__.get("dem", None)
        if isinstance(dem, str):
            dic.update({"dem": dem})
        elif isinstance(dem, pdem.Dem):
            dic.update({"dem": dem.Dataset.elevation.attrs})

        meteo = self.__dict__.get("meteo", None)
        if isinstance(meteo, str):
            dic.update({"meteo": meteo})
        elif isinstance(meteo, pmeteo.Meteo):
            dic.update({"meteo": meteo.Dataset.attrs})

        dic["version"] = pyposeidon.__version__

        for attr, value in dic.items():
            if isinstance(value, datetime.datetime):
                dic[attr] = dic[attr].isoformat()
            if isinstance(value, pd.Timedelta):
                dic[attr] = dic[attr].isoformat()
            if isinstance(value, pd.DataFrame):
                dic[attr] = dic[attr].to_dict()
        json.dump(dic, open(path + self.tag + "_model.json", "w"), default=myconverter)

    def output(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./d3d/")
        slevel = get_value(self, kwargs, "slevel", 0.0)
        flag = get_value(self, kwargs, "update", [])

        nj, ni = self.mesh.Dataset.lons.shape

        if not os.path.exists(path):
            os.makedirs(path)

        # save mdf
        self.mdf.to_csv(path + self.tag + ".mdf", sep="=")

        # save mesh file
        if flag:
            if ("all" in flag) | ("mesh" in flag):
                # save mesh
                self.mesh.to_file(filename=path + self.tag + ".grd")
            else:
                logger.info("skipping mesh file ..\n")
        else:
            self.mesh.to_file(filename=path + self.tag + ".grd")

        # save bathymetry file
        self.to_dep(self.dem.Dataset, rpath=path, tag=self.tag, update=flag)

        # save meteo
        if self.atm:
            try:
                self.to_force(self.meteo.Dataset, vars=["msl", "u10", "v10"], rpath=path, **kwargs)
            except AttributeError as e:
                logger.warning("no meteo data available.. no update..\n")
                pass

        # save obs file
        self.to_obs(self, **kwargs)

        # save enc file
        if flag:

            if ("all" in flag) | ("model" in flag):
                # save enc
                # write enc out
                with open(path + self.tag + ".enc", "w") as f:
                    f.write("{:>5}{:>5}\n".format(ni + 1, 1))  # add one like ddb
                    f.write("{:>5}{:>5}\n".format(ni + 1, nj + 1))
                    f.write("{:>5}{:>5}\n".format(1, nj + 1))
                    f.write("{:>5}{:>5}\n".format(1, 1))
                    f.write("{:>5}{:>5}\n".format(ni + 1, 1))

        else:

            # write enc out
            with open(path + self.tag + ".enc", "w") as f:
                f.write("{:>5}{:>5}\n".format(ni + 1, 1))  # add one like ddb
                f.write("{:>5}{:>5}\n".format(ni + 1, nj + 1))
                f.write("{:>5}{:>5}\n".format(1, nj + 1))
                f.write("{:>5}{:>5}\n".format(1, 1))
                f.write("{:>5}{:>5}\n".format(ni + 1, 1))

        calc_dir = get_value(self, kwargs, "rpath", "./d3d/")

        try:
            bin_path = os.environ["D3D"]
        except:
            bin_path = get_value(self, kwargs, "epath", None)

        try:
            lib_path = os.environ["LD3D"]
        except:
            lib_path = get_value(self, kwargs, "lpath", None)

        if bin_path is None:
            # ---------------------------------------------------------------------
            logger.warning("D3D executable path (epath) not given\n")
            # ---------------------------------------------------------------------

        if lib_path is None:
            # ---------------------------------------------------------------------
            logger.warning("D3D libraries path (lpath) not given\n")
            # ---------------------------------------------------------------------

        ncores = get_value(self, kwargs, "ncores", NCORES)

        if not os.path.exists(calc_dir + self.tag + "_hydro.xml"):

            # edit and save config file
            copy2(DATA_PATH + "config_d_hydro.xml", calc_dir + self.tag + "_hydro.xml")

        xml = md.parse(calc_dir + self.tag + "_hydro.xml")

        xml.getElementsByTagName("mdfFile")[0].firstChild.replaceWholeText(self.tag + ".mdf")

        with open(calc_dir + self.tag + "_hydro.xml", "w") as f:
            xml.writexml(f)

        if not os.path.exists(calc_dir + "run_flow2d3d.sh"):

            copy2(DATA_PATH + "run_flow2d3d.sh", calc_dir + "run_flow2d3d.sh")

            # make the script executable
            execf = calc_dir + "run_flow2d3d.sh"
            mode = os.stat(execf).st_mode
            mode |= (mode & 0o444) >> 2  # copy R bits to X
            os.chmod(execf, mode)

        # ---------------------------------------------------------------------
        logger.info("output done\n")
        # ---------------------------------------------------------------------

    @staticmethod
    def vpoints(grid, obs_points, bat, **kwargs):

        idx = []
        jdx = []
        for m in range(obs_points.shape[0]):
            lat, lon = obs_points.loc[m, ["lat", "lon"]]
            nearest = grid.sel(x=[lon], y=[lat], method="nearest")
            j = np.abs(grid.x.values - nearest.x.values).argmin()
            i = np.abs(grid.y.values - nearest.y.values).argmin()
            if bat[i, j]:
                idx.append(i)
                jdx.append(j)
            else:
                bnear = bat[i - 5 : i + 6, j - 5 : j + 6]  # near by grid nodes

                rlon = grid.lons[i - 5 : i + 6, j - 5 : j + 6] - lon
                rlat = grid.lats[i - 5 : i + 6, j - 5 : j + 6] - lat
                rad = np.sqrt(rlon ** 2 + rlat ** 2)  # radial distance from the obs point

                rmask = rad.values[bnear.mask == False]  # mask the distance array with the valid mask from dem

                rmask.sort()  # sort to start close and move further away
                if rmask.size > 0:

                    for r in rmask:  # Find the closest valid point
                        [[k, l]] = np.argwhere(rad.values == r)
                        if bnear[k - 1 : k + 1, l - 1 : l + 1].mask.sum() == 0:
                            break  # The point is valid point

                    xv = rad[k, l].x.values  # lat, lon of valid point
                    yv = rad[k, l].y.values

                    # final i,j
                    j = np.abs(grid.x.values - xv).argmin()
                    i = np.abs(grid.y.values - yv).argmin()

                    idx.append(i)
                    jdx.append(j)

                else:

                    idx.append(np.nan)
                    jdx.append(np.nan)

        return idx, jdx

    def execute(self, **kwargs):

        self.create(**kwargs)
        self.output(**kwargs)
        self.save(**kwargs)
        self.run(**kwargs)

    def read_folder(self, rfolder, **kwargs):

        gfile = glob.glob(rfolder + "/*.grd")  # Grid
        dfile = glob.glob(rfolder + "/*.dep")  # bathymetry
        u = glob.glob(rfolder + "/*.amu")  # meteo
        v = glob.glob(rfolder + "/*.amv")
        p = glob.glob(rfolder + "/*.amp")

        # config
        self.mdf = pd.read_csv(d[0], sep="=")
        self.mdf = self.mdf.set_index(self.mdf.columns[0])  # set index
        # mesh
        self.mesh = pmesh.set("r2d", mesh_file=gfile[0])
        # bath
        self.dem.Dataset = d3d.from_dep(dfile[0])
        # meteo
        mf = []
        mf.append(d3d.from_force(u[0], "u10"))
        mf.append(d3d.from_force(v[0], "v10"))
        mf.append(d3d.from_force(p[0], "msl"))
        self.meteo.Dataset = xr.merge(mf)

        # ---------------------------------------------------------------------
        logger.exception("folder incomplete. Abort\n")
        sys.exit(1)
        # ---------------------------------------------------------------------

    def get_data(self, **kwargs):

        dic = self.__dict__

        dic.update(kwargs)

        self.data = data.get_output(**dic)
