"""
Schism model of pyposeidon. It controls the creation, output & execution of a complete simulation based on schism

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
import pkg_resources
import json
from collections import OrderedDict
import pandas as pd
import geopandas as gp
import glob
from shutil import copyfile
import xarray as xr
import cartopy.feature as cf
import geopandas as gp
import f90nml
import errno
from tqdm import tqdm

# local modules
import pyposeidon
import pyposeidon.grid as pgrid
import pyposeidon.meteo as pmeteo
import pyposeidon.dem as pdem
from pyposeidon.utils.get_value import get_value
from pyposeidon.utils.converter import myconverter
from pyposeidon.utils.vals import obs
from pyposeidon.utils.cpoint import closest_node
from pyposeidon.utils.hfun import hfun
from pyposeidon.utils.unml import unml
from pyposeidon.utils.data import data

import logging

logger = logging.getLogger("pyposeidon")

import multiprocessing

NCORES = max(1, multiprocessing.cpu_count() - 1)

# retrieve the module path
# DATA_PATH = pkg_resources.resource_filename('pyposeidon', 'misc')
DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/misc/"
TEST_DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/tests/data/"

# add conda path to PATH
cpath = pyposeidon.__path__[0].split("/lib/")[0]
os.environ["PATH"] += os.pathsep + cpath + "/bin"


class schism:
    def __init__(self, **kwargs):
        """
        Create a Schism solver

        !!! danger ""
            Due to a limitation of the Library rendering the docstrings, all arguments are marked
            as `required`, nevertheless they are all `Optional`.

        Args:
            rfolder (str): The path to a directory containing the results of a Model solved by Schism.
            global_grid (str):
            geometry (Union[dict, str]): A dictionary containing the the bounding box of the region
                you want to solve or a shapefile.
            load_grid (bool): Flag indicating whether to load the grid or not. Defaults to `False`.
            load_meteo (bool): Flag indicating whether to load the meteo data or not. Defauls to
                `False`.
            coastlines (Union[str, GeoDataFrame]): A `GeoDataFrame` or the path to a shapefile which
                describes the coastlines.
            coast_resolution (str): If no coastlines have been provided, then Cartopy's
                `NaturalEarthFeatures` are being used.  The resolution of the coastlines can be
                controlled with this parameter.  Valid values are: `{"l", "i", "h"}` corresponding to
                resolutions of `{110, 50, 10}` meters respectively. Defaults to `i`.
            tag: (str): The solver's "tag". Defaults to `"schism"`. **Do we need this?**
            tide: (str): Flag indicating whether to load "tide". Defaults to `False`.
            atm (bool): The solver's atm. Defaults to `True`.
            monitor (bool: The solver's monitor. Defaults to `False`.
            epath (str): The path to the schism executable. If the `SCHISM` env variable has been
                set, then it overrides the value passed as the parameter.
            start_date (str): The date from which the analysis should start. It should be a string parseable
                by `pd.to_datetime()`.
            end_date (str): The date at which the analysis should end. It should be a string parseable by
                `pd.to_datetime()`.
            time_frame (str): The duration of the analysis. It should be a string parseable by
                `pd.to_datetime()`.
        """

        rfolder = kwargs.get("rfolder", None)
        if rfolder:
            self.read_folder(**kwargs)

        self.global_grid = kwargs.get("global_grid", False)
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

                self.lon_min, self.lat_min, self.lon_max, self.lat_max = geo.total_bounds

        # coastlines
        coastlines = kwargs.get("coastlines", None)

        if coastlines is None:
            cr = kwargs.get("coast_resolution", "i")

            # world polygons - user input
            coast = cf.NaturalEarthFeature(
                category="physical", name="land", scale="{}m".format({"l": 110, "i": 50, "h": 10}[cr])
            )

            self.coastlines = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])

        else:

            try:
                coast = gp.GeoDataFrame.from_file(coastlines)
            except:
                coast = gp.GeoDataFrame(coastlines)

            self.coastlines = coast

        if not hasattr(self, "start_date"):
            start_date = kwargs.get("start_date", None)
            self.start_date = pd.to_datetime(start_date)

        if not hasattr(self, "end_date"):

            if "time_frame" in kwargs:
                time_frame = kwargs.get("time_frame", None)
                self.end_date = self.start_date + pd.to_timedelta(time_frame)
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

        self.tag = kwargs.get("tag", "schism")
        self.tide = kwargs.get("tide", False)
        self.atm = kwargs.get("atm", True)
        self.monitor = kwargs.get("monitor", False)

        try:
            self.epath = os.environ["SCHISM"]
        except:
            self.epath = kwargs.get("epath", None)

        self.solver = self.__class__.__name__

        for attr, value in kwargs.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

        setattr(self, "misc", {})

    # ============================================================================================
    # CONFIG
    # ============================================================================================

    def config(self, config_file=None, output=False, **kwargs):

        dic = get_value(self, kwargs, "parameters", None)
        #        param_file = get_value(self,kwargs,'config_file',None)

        if config_file:
            # ---------------------------------------------------------------------
            logger.info("reading parameter file {}\n".format(config_file))
            # ---------------------------------------------------------------------
        else:
            # ---------------------------------------------------------------------
            logger.info("using default parameter file ...\n")
            # ---------------------------------------------------------------------

            config_file = DATA_PATH + "param.nml"

        params = f90nml.read(config_file)

        # update key values

        patch = {
            "start_year": self.start_date.year,
            "start_month": self.start_date.month,
            "start_day": self.start_date.day,
            "start_hour": self.start_date.hour,
        }

        params = unml(params, patch)
        # update
        if dic:
            params = unml(params, dic)

        # test rnday
        if float(params["CORE"]["rnday"]) * 24 * 3600 > (self.end_date - self.start_date).total_seconds():
            # ---------------------------------------------------------------------
            logger.warning("rnday larger than simulation range\n")
            logger.warning(
                "rnday={} while simulation time is {}\n".format(
                    params["core"]["rnday"], (self.end_date - self.start_date).total_seconds() / (3600 * 24.0)
                )
            )
            # ---------------------------------------------------------------------

        self.params = params

        if output:
            # save params
            # ---------------------------------------------------------------------
            logger.info("output param.nml file ...\n")
            # ---------------------------------------------------------------------

            path = get_value(self, kwargs, "rpath", "./schism/")
            self.params.write(path + "param.nml", force=True)

    # ============================================================================================
    # METEO
    # ============================================================================================
    def force(self, **kwargs):

        meteo_source = get_value(self, kwargs, "meteo_source", None)

        kwargs.update({"meteo_source": meteo_source})

        flag = get_value(self, kwargs, "update", ["all"])

        z = {**self.__dict__, **kwargs}  # merge self and possible kwargs

        # check if files exist
        if flag:
            if ("meteo" in flag) | ("all" in flag):
                self.meteo = pmeteo.meteo(**z)
            else:
                logger.info("skipping meteo ..\n")
        else:
            self.meteo = pmeteo.meteo(**z)

        if hasattr(self, "meteo"):
            # add 1 hour for Schism issue with end time
            ap = self.meteo.Dataset.isel(time=-1)
            ap["time"] = ap.time.values + pd.to_timedelta("1H")

            self.meteo.Dataset = xr.concat([self.meteo.Dataset, ap], dim="time")

    @staticmethod
    def to_force(ar0, **kwargs):

        logger.info("writing meteo files ..\n")

        path = kwargs.get("rpath", "./schism/")

        [p, u, v] = kwargs.get("vars", "[None,None,None]")

        ar = ar0.sortby("latitude", ascending=True)

        xx, yy = np.meshgrid(ar.longitude.data, ar.latitude.data)

        zero = np.zeros(ar[p].data.shape)

        date = kwargs.get("date", ar.time[0].data)

        udate = pd.to_datetime(date).strftime("%Y-%m-%d")

        bdate = pd.to_datetime(date).strftime("%Y %m %d %H").split(" ")

        tlist = (ar.time.data - pd.to_datetime([udate]).values).astype("timedelta64[s]") / 3600.0

        tlist = tlist.astype(float) / 24.0

        bdate = [int(q) for q in bdate[:3]] + [0]

        sout = xr.Dataset(
            {
                "prmsl": (["time", "nx_grid", "ny_grid"], ar[p].data),
                "uwind": (["time", "nx_grid", "ny_grid"], ar[u].data),
                "vwind": (["time", "nx_grid", "ny_grid"], ar[v].data),
                "spfh": (["time", "nx_grid", "ny_grid"], zero),
                "stmp": (["time", "nx_grid", "ny_grid"], zero),
                "lon": (["nx_grid", "ny_grid"], xx),
                "lat": (["nx_grid", "ny_grid"], yy),
            },
            coords={"time": tlist},
        )

        sout.attrs = {"description": "Schism meteo data", "history": "pyposeidon", "source": "netCDF4 python module"}

        sout.time.attrs = {"long_name": "Time", "standard_name": "time", "base_date": bdate, "units": udate}

        sout.lat.attrs = {"units": "degrees_north", "long_name": "Latitude", "standard_name": "latitude"}

        sout.prmsl.attrs = {
            "units": "Pa",
            "long_name": "Pressure reduced to MSL",
            "standard_name": "air_pressure_at_sea_level",
        }

        sout.uwind.attrs = {
            "units": "m/s",
            "long_name": "Surface Eastward Air Velocity",
            "standard_name": "eastward_wind",
        }

        sout.vwind.attrs = {
            "units": "m/s",
            "long_name": "Surface Northward Air Velocity",
            "standard_name": "northward_wind",
        }

        sout.spfh.attrs = {
            "units": "1",
            "long_name": "Surface Specific Humidity (2m AGL)",
            "standard_name": "specific_humidity",
        }

        sout.stmp.attrs = {
            "units": "degrees",
            "long_name": "Surface Temperature",
            "standard_name": "surface temperature",
        }

        # check if folder sflux exists
        if not os.path.exists(path + "sflux"):
            os.makedirs(path + "sflux")

        m_index = kwargs.get("m_index", 1)

        filename = kwargs.get("filename", "sflux/sflux_air_{}.0001.nc".format(m_index))

        sout.to_netcdf(path + filename)

    # ============================================================================================
    # DEM
    # ============================================================================================

    def bath(self, **kwargs):
        #       z = self.__dict__.copy()

        kwargs["grid_x"] = self.grid.Dataset.SCHISM_hgrid_node_x.values
        kwargs["grid_y"] = self.grid.Dataset.SCHISM_hgrid_node_y.values

        dpath = get_value(self, kwargs, "dem_source", None)

        kwargs.update({"dem_source": dpath})

        flag = get_value(self, kwargs, "update", ["all"])
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
                self.dem = pdem.dem(**kwargs)
            else:
                logger.info("dem from grid file\n")

    # ============================================================================================
    # EXECUTION
    # ============================================================================================
    def create(self, **kwargs):

        if not kwargs:
            kwargs = self.__dict__.copy()

        # Grid

        self.grid = pgrid.grid(type="tri2d", **kwargs)

        # set lat/lon from file
        if hasattr(self, "grid_file"):
            kwargs.update({"lon_min": self.grid.Dataset.SCHISM_hgrid_node_x.values.min()})
            kwargs.update({"lon_max": self.grid.Dataset.SCHISM_hgrid_node_x.values.max()})
            kwargs.update({"lat_min": self.grid.Dataset.SCHISM_hgrid_node_y.values.min()})
            kwargs.update({"lat_max": self.grid.Dataset.SCHISM_hgrid_node_y.values.max()})

            self.lon_min = self.grid.Dataset.SCHISM_hgrid_node_x.values.min()
            self.lon_max = self.grid.Dataset.SCHISM_hgrid_node_x.values.max()
            self.lat_min = self.grid.Dataset.SCHISM_hgrid_node_y.values.min()
            self.lat_max = self.grid.Dataset.SCHISM_hgrid_node_y.values.max()

        # get bathymetry
        self.bath(**kwargs)

        # get boundaries
        # self.bc()

        # get meteo
        if self.atm:
            self.force(**kwargs)

        # get tide
        if self.tide:
            self.tidebc()

        self.config(**kwargs)

    def output(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")
        flag = get_value(self, kwargs, "update", ["all"])
        split_by = get_value(self, kwargs, "meteo_split_by", None)

        if not os.path.exists(path):
            os.makedirs(path)

        # save sflux_inputs.txt
        if not os.path.exists(path + "sflux"):
            os.makedirs(path + "sflux")

        with open(path + "sflux/sflux_inputs.txt", "w") as f:
            f.write("&sflux_inputs\n")
            f.write("/ \n\n")

        # save bctides.in
        bs = self.grid.Dataset[["node", "id", "type"]].to_dataframe()
        # open boundaries
        number_of_open_boundaries = bs.id.max()
        number_of_open_boundaries_nodes = bs.loc[bs.id > 0].shape[0]

        with open(path + "bctides.in", "w") as f:
            f.write("Header\n")
            f.write("{} {}\n".format(0, 40.0))  #  ntip tip_dp
            f.write("{}\n".format(0))  # nbfr
            f.write("{}\n".format(number_of_open_boundaries))  # number of open boundaries
            for i in range(1, number_of_open_boundaries + 1):
                nnodes = bs.loc[bs.id == i, "node"].shape[0]
                f.write(
                    "{} {} {} {} {}\n".format(nnodes, 2, 0, 0, 0)
                )  # number of nodes on the open boundary segment j (corresponding to hgrid.gr3), B.C. flags for elevation, velocity, temperature, and salinity
                f.write("{}\n".format(0))  # ethconst !constant elevation value for this segment

        # save vgrid.in
        with open(path + "vgrid.in", "w") as f:
            f.write("{}\n".format(2))  # ivcor (1: LSC2; 2: SZ)
            f.write(
                "{} {} {}\n".format(2, 1, 1.0e6)
            )  # nvrt(=Nz); kz (# of Z-levels); hs (transition depth between S and Z)
            f.write("Z levels\n")  # Z levels !Z-levels in the lower portion
            f.write(
                "{} {}\n".format(1, -1.0e6)
            )  #!level index, z-coordinates, z-coordinate of the last Z-level must match -hs
            f.write("S levels\n")  # S-levels below
            f.write(
                "{} {} {}\n".format(40.0, 1.0, 1.0e-4)
            )  # constants used in S-transformation: h_c, theta_b, theta_f
            f.write("{} {}\n".format(1, -1.0))  # first S-level (sigma-coordinate must be -1)
            f.write("{} {}\n".format(2, 0.0))  # levels index, sigma-coordinate, last sigma-coordinate must be 0

        # save params.in

        self.params.write(path + "param.nml", force=True)

        # save hgrid.gr3
        try:

            try:

                bat = -self.dem.Dataset.fval.values.astype(float)  # minus for the hydro run

            except:

                bat = -self.dem.Dataset.ival.values.astype(float)  # minus for the hydro run

            self.grid.Dataset.depth.loc[: bat.size] = bat

            self.grid.to_file(filename=path + "hgrid.gr3")
            copyfile(path + "hgrid.gr3", path + "hgrid.ll")

            logger.info("updating bathymetry ..\n")

        except AttributeError as e:

            logger.info("Keeping bathymetry from hgrid.gr3 ..\n")

            copyfile(self.grid_file, path + "hgrid.gr3")  # copy original grid file
            copyfile(path + "hgrid.gr3", path + "hgrid.ll")

        # manning file
        manfile = path + "manning.gr3"

        if hasattr(self, "manning_file"):
            copyfile(self.manning_file, manfile)  # copy original manning file
            if self.manning_file == manfile:
                logger.info("Keeping manning file ..\n")

        manning = get_value(self, kwargs, "manning", 0.12)
        nn = self.grid.Dataset.nSCHISM_hgrid_node.size
        n3e = self.grid.Dataset.nSCHISM_hgrid_face.size

        with open(manfile, "w") as f:
            f.write("\t 0 \n")
            f.write("\t {} {}\n".format(n3e, nn))

        df = self.grid.Dataset[["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]].to_dataframe()

        df["man"] = manning

        df.index = np.arange(1, len(df) + 1)

        df.to_csv(
            manfile,
            index=True,
            sep="\t",
            header=None,
            mode="a",
            float_format="%.10f",
            columns=["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "man"],
        )

        logger.info("Manning file created..\n")

        # windrot_geo2proj

        windfile = path + "windrot_geo2proj.gr3"

        if hasattr(self, "windrot_file"):
            copyfile(self.windrot_file, windfile)  # copy original grid file
            if self.windrot_file != windfile:
                logger.info("Keeping windrot_geo2proj file ..\n")

        windrot = get_value(self, kwargs, "windrot", 0.00001)

        with open(windfile, "w") as f:
            f.write("\t 0 \n")
            f.write("\t {} {}\n".format(n3e, nn))

        df = self.grid.Dataset[["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]].to_dataframe()

        df["windrot"] = windrot

        df.index = np.arange(1, len(df) + 1)

        df.to_csv(
            windfile,
            index=True,
            sep="\t",
            header=None,
            mode="a",
            float_format="%.10f",
            columns=["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "windrot"],
        )

        logger.info("Windrot_geo2proj file created..\n")

        # save meteo
        m_index = get_value(self, kwargs, "m_index", 1)

        if hasattr(self, "atm"):
            try:
                if split_by:
                    times, datasets = zip(*self.meteo.Dataset.groupby("time.{}".format(split_by)))
                    mpaths = ["sflux/sflux_air_{}.{:04d}.nc".format(m_index, t + 1) for t in np.arange(len(times))]
                    for das, mpath in list(zip(datasets, mpaths)):
                        self.to_force(
                            das, vars=["msl", "u10", "v10"], rpath=path, filename=mpath, date=self.date, **kwargs
                        )
                else:
                    self.to_force(self.meteo.Dataset, vars=["msl", "u10", "v10"], rpath=path, **kwargs)
            except AttributeError as e:
                logger.warning("no meteo data available.. no update..\n")
                pass

        calc_dir = get_value(self, kwargs, "rpath", "./schism/")

        try:
            bin_path = os.environ["SCHISM"]
        except:
            bin_path = get_value(self, kwargs, "epath", None)

        if bin_path is None:
            # ------------------------------------------------------------------------------
            logger.warning("Schism executable path (epath) not given -> using default \n")
            # ------------------------------------------------------------------------------
            bin_path = "schism"

        ncores = get_value(self, kwargs, "ncores", NCORES)

        with open(calc_dir + "launchSchism.sh", "w") as f:
            f.write("exec={}\n".format(bin_path))
            f.write("mkdir outputs\n")
            f.write("mpirun -N {} $exec\n".format(ncores))

        # make the script executable
        execf = calc_dir + "launchSchism.sh"
        mode = os.stat(execf).st_mode
        mode |= (mode & 0o444) >> 2  # copy R bits to X
        os.chmod(execf, mode)

        # ---------------------------------------------------------------------
        logger.info("output done\n")
        # ---------------------------------------------------------------------

    def run(self, **kwargs):

        calc_dir = get_value(self, kwargs, "rpath", "./schism/")

        ncores = get_value(self, kwargs, "ncores", NCORES)

        # ---------------------------------------------------------------------
        logger.info("executing model\n")
        # ---------------------------------------------------------------------

        # note that cwd is the folder where the executable is
        ex = subprocess.Popen(
            args=["./launchSchism.sh"], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
        )  # , bufsize=1)

        with open(calc_dir + "err.log", "w") as f:
            for line in iter(ex.stderr.readline, b""):
                f.write(line.decode(sys.stdout.encoding))
                logger.info(line.decode(sys.stdout.encoding))
        ex.stderr.close()

        with open(calc_dir + "run.log", "w") as f:
            for line in iter(ex.stdout.readline, b""):
                f.write(line.decode(sys.stdout.encoding))
                logger.info(line.decode(sys.stdout.encoding))
        ex.stdout.close()

        # ---------------------------------------------------------------------
        logger.info("FINISHED\n")
        # ---------------------------------------------------------------------

    def save(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")

        lista = [key for key, value in self.__dict__.items() if key not in ["meteo", "dem", "grid"]]
        dic = {k: self.__dict__.get(k, None) for k in lista}

        grid = self.__dict__.get("grid", None)
        if isinstance(grid, str):
            dic.update({"grid": grid})
        else:
            dic.update({"grid": grid.__class__.__name__})

        dem = self.__dict__.get("dem", None)
        if isinstance(dem, str):
            dic.update({"dem": dem})
        elif isinstance(dem, pdem.dem):
            dic.update({"dem": dem.Dataset.elevation.attrs})

        meteo = self.__dict__.get("meteo", None)
        if isinstance(meteo, str):
            dic.update({"meteo": meteo})
        elif isinstance(meteo, pmeteo.meteo):
            try:
                dic.update({"meteo": [meteo.Dataset.attrs]})
            except:
                dic.update({"meteo": [x.attrs for x in meteo.Dataset]})

        coast = self.__dict__.get("coast_resolution", None)
        coastline = self.__dict__.get("coastlines", None)
        if isinstance(coast, str):
            dic.update({"coastlines": None})
        elif isinstance(coastline, gp.GeoDataFrame):  # TODO
            dic.update({"coastlines": None})

        dic["version"] = pyposeidon.__version__

        for attr, value in dic.items():
            if isinstance(value, datetime.datetime):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.Timedelta):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.DataFrame):
                dic[attr] = value.to_dict()

        json.dump(dic, open(path + self.tag + "_model.json", "w"), default=myconverter)

    def execute(self, **kwargs):

        flag = get_value(self, kwargs, "update", ["all"])
        if flag:
            self.create(**kwargs)
            self.output(**kwargs)
            if self.monitor:
                self.set_obs()
            self.save(**kwargs)
            self.run(**kwargs)

        else:

            self.save(**kwargs)

            calc_dir = get_value(self, kwargs, "rpath", "./schism/")

            try:
                bin_path = os.environ["SCHISM"]
            except:
                bin_path = get_value(self, kwargs, "epath", None)

            if bin_path is None:
                # ------------------------------------------------------------------------------
                logger.warning("Schism executable path (epath) not given -> using default \n")
                # ------------------------------------------------------------------------------
                bin_path = "schism"

            ncores = get_value(self, kwargs, "ncores", NCORES)

            with open(calc_dir + "launchSchism.sh", "w") as f:
                f.write("exec={}\n".format(bin_path))
                f.write("mkdir outputs\n")
                f.write("mpirun -N {} $exec\n".format(ncores))

            # make the script executable
            execf = calc_dir + "launchSchism.sh"
            mode = os.stat(execf).st_mode
            mode |= (mode & 0o444) >> 2  # copy R bits to X
            os.chmod(execf, mode)

            self.run(**kwargs)

    def read_folder(self, rfolder, **kwargs):

        self.rpath = rfolder
        s = glob.glob(rfolder + "/param.nml")
        mfiles1 = glob.glob(rfolder + "/sflux/*_1*.nc")
        mfiles1.sort()
        # check for 2nd meteo
        mfiles2 = glob.glob(rfolder + "/sflux/*_2*.nc")
        mfiles2.sort()

        mfiles = {"1": mfiles1, "2": mfiles2}

        mfiles = {k: v for k, v in mfiles.items() if v}  # remove empty keys, e.g. no mfiles2

        hfile = rfolder + "/hgrid.gr3"  # Grid
        self.params = f90nml.read(s[0])

        mykeys = ["start_year", "start_month", "start_day"]
        sdate = [self.params["opt"][x] for x in mykeys]  # extract date info from param.nml
        sd = "-".join(str(item) for item in sdate)  # join in str
        sd = pd.to_datetime(sd)  # convert to datestamp

        sh = [self.params["opt"][x] for x in ["start_hour", "utc_start"]]  # get time attrs
        sd = sd + pd.to_timedelta("{}H".format(sh[0] + sh[1]))  # complete start_date

        ed = sd + pd.to_timedelta("{}D".format(self.params["core"]["rnday"]))  # compute end date based on rnday

        self.start_date = sd  # set attrs
        self.end_date = ed

        load_grid = get_value(self, kwargs, "load_grid", False)

        if load_grid:
            try:
                self.grid = pgrid.grid(type="tri2d", grid_file=hfile)
            except:
                logger.warning("Loading grid failed")
                pass
        else:
            logger.warning("No grid loaded")

        load_meteo = get_value(self, kwargs, "load_meteo", False)

        # meteo
        if load_meteo is True:

            try:
                pm = []
                for key, val in mfiles.items():
                    msource = xr.open_mfdataset(val)
                    pm.append(msource)
                if len(pm) == 1:
                    pm = pm[0]
                self.meteo = pmeteo.meteo(meteo_source=pm, meteo_engine="passthrough")

            except:
                pm = []
                for key, val in mfiles.items():
                    ma = []

                    for ifile in mfiles:
                        g = xr.open_dataset(ifile)
                        ts = "-".join(g.time.attrs["base_date"].astype(str)[:3])
                        time_r = pd.to_datetime(ts)
                        try:
                            times = time_r + pd.to_timedelta(g.time.values, unit="D").round("H")
                            g = g.assign_coords({"time": times})
                            ma.append(g)
                        except:
                            logger.warning("Loading meteo failed")
                            break
                    if ma:
                        msource = xr.merge(ma)
                        pm.append(msource)

                if len(pm) == 1:
                    pm = pm[0]

                self.meteo = pmeteo.meteo(meteo_source=pm, meteo_engine="passthrough")

        else:
            logger.warning("No meteo loaded")

    def global2local(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")

        # Read the global node index distribution to the cores
        gfiles = glob.glob(path + "outputs/local_to_global_*")
        gfiles.sort()

        # create a dict from filenames to identify parts in the dataframes below
        keys = []
        for name in gfiles:
            keys.append("core{}".format(name.split("/")[-1].split("_")[-1]))

        # Parsing the files

        # We read from the first file the header (it is the same for all)
        with open(gfiles[0], "r") as f:
            header = pd.read_csv(
                f,
                header=None,
                nrows=1,
                delim_whitespace=True,
                names=[
                    "ns_global",
                    "ne_global",
                    "np_global",
                    "nvrt",
                    "nproc",
                    "ntracers",
                    "T",
                    "S",
                    "GEN",
                    "AGE",
                    "SED3D",
                    "EcoSim",
                    "ICM",
                    "CoSINE",
                    "Feco",
                    "TIMOR",
                    "FABM",
                    "DVD",
                ],
            )

        # get the number of elems from all files
        nels = []
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                ne = pd.read_csv(f, skiprows=2, header=None, nrows=1)
                nels.append(ne.values.flatten()[0].astype(int))

        # read and add them to pandas DataFrame
        frames = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                frames[i] = pd.read_csv(
                    f, skiprows=3, header=None, nrows=nels[i], names=["local", "global_n"], delim_whitespace=True
                )

        elems = pd.concat(frames, keys=keys)

        # get the number of nodes from all files
        nq = []
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                nn = pd.read_csv(f, skiprows=nels[i] + 3, header=None, nrows=1)
                nq.append(nn.values.flatten()[0].astype(int))

        # read and add them to pandas DataFrame
        nframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                nframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + 4,
                    header=None,
                    nrows=nq[i],
                    names=["local", "global_n"],
                    delim_whitespace=True,
                )

        nodes = pd.concat(nframes, keys=keys)

        # get the number of edges
        nw = []
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                nb = pd.read_csv(f, skiprows=nels[i] + nq[i] + 4, header=None, nrows=1)
                nw.append(nb.values.flatten()[0].astype(int))

        # read and add them to pandas DataFrame
        wframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                wframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + nq[i] + 5,
                    header=None,
                    nrows=nw[i],
                    names=["local", "global_n"],
                    delim_whitespace=True,
                )

        re = pd.concat(wframes, keys=keys)

        # read secondary headers
        with open(gfiles[0], "r") as f:
            h0 = pd.read_csv(
                f,
                skiprows=nels[0] + nq[0] + nw[0] + 6,
                header=None,
                nrows=1,
                delim_whitespace=True,
                names=["start_year", "start_month", "start_day", "start_hour", "utc_start"],
            )
        with open(gfiles[0], "r") as f:
            h1 = pd.read_csv(
                f,
                skiprows=nels[0] + nq[0] + nw[0] + 7,
                header=None,
                nrows=1,
                delim_whitespace=True,
                names=["nrec", "dtout", "nspool", "nvrt", "kz", "h0", "h_s", "h_c", "theta_b", "theta_f", "ics"],
            )

        ztots = ["ztot_" + str(i) for i in range(1, h1.loc[:, "kz"].values[0] - 1)]
        sigmas = ["sigma_" + str(i) for i in range(h1.loc[:, "nvrt"].values[0] - h1.loc[:, "kz"].values[0] + 1)]

        # read secondary header
        with open(gfiles[0], "r") as f:
            h2 = pd.read_csv(
                f,
                skiprows=nels[0] + nq[0] + nw[0] + 8,
                header=None,
                nrows=1,
                delim_whitespace=True,
                names=ztots + sigmas,
            )

        # combine headers
        self.misc.update({"header": pd.concat([h0, h1, h2], axis=1)})

        # read lat/lon from all files
        gframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                gframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + nq[i] + nw[i] + 10,
                    header=None,
                    nrows=nq[i],
                    delim_whitespace=True,
                    names=["lon", "lat", "depth", "kbp00"],
                )

        grid = pd.concat(gframes, keys=keys)

        # Droping duplicates
        drops = nodes.reset_index()[nodes.reset_index().duplicated("global_n")].index.to_list()
        cnodes = nodes.global_n.drop_duplicates()  # drop duplicate global nodes and store the values to an array

        grid = grid.reset_index().drop(drops)  # Use the mask from nodes to match grid
        grid = grid.set_index(["level_0", "level_1"])

        grid.index = grid.index.droplevel()  # drop multi-index
        grid = grid.reset_index(drop=True)  # reset index
        grid.index = cnodes.values - 1  # reindex based on the global index, -1 for the python convention
        grd = grid.sort_index()  # sort with the new index (that is the global_n)
        self.misc.update({"grd": grd.reset_index(drop=True)})  # reindex for final version

        # Read tessalation
        eframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                eframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + nq[i] + nw[i] + nq[i] + 10,
                    header=None,
                    nrows=nels[i],
                    delim_whitespace=True,
                    names=["type", "a", "b", "c", "d"],
                )

        tri = pd.concat(eframes, keys=keys)

        # Replace local index with the global one
        for key in keys:
            nod = nodes.loc[key].copy()  # make a copy of the the core
            nod = nod.set_index(
                "local"
            )  # reset the index using the local values (in essense set the index to start from 1...)
            tri.loc[key, "ga"] = nod.reindex(tri.loc[key, "a"].values).values
            tri.loc[key, "gb"] = nod.reindex(tri.loc[key, "b"].values).values
            tri.loc[key, "gc"] = nod.reindex(tri.loc[key, "c"].values).values
            tri.loc[key, "gd"] = nod.reindex(tri.loc[key, "d"].values).values

        tri = tri.drop(["a", "b", "c", "d"], axis=1)  # drop local references

        # sort
        gt3 = tri.loc[:, ["type", "ga", "gb", "gc", "gd"]].copy()  # make a copy
        gt3.index = gt3.index.droplevel()  # drop multi-index
        gt3 = gt3.reset_index(drop=True)
        # Now we need to put them in order based on the global index in elems
        gt3.index = elems.global_n.values  # we set the index equal to the global_n column
        gt3 = gt3.sort_index()  # sort them

        # add nan column in place of the fourth node. NOTE:  This needs to be tested for quadrilaterals
        #        gt3['gd']=np.nan

        gt3 = gt3.reset_index()  # reset to add more columns without problems

        ## Add mean x, y of the elememts. To be used in the output
        gt3["x1"] = grd.loc[gt3["ga"].values - 1, "lon"].values  # lon of the index, -1 for python convention
        gt3["y1"] = grd.loc[gt3["ga"].values - 1, "lat"].values  # lat of the index
        gt3["x2"] = grd.loc[gt3["gb"].values - 1, "lon"].values
        gt3["y2"] = grd.loc[gt3["gb"].values - 1, "lat"].values
        gt3["x3"] = grd.loc[gt3["gc"].values - 1, "lon"].values
        gt3["y3"] = grd.loc[gt3["gc"].values - 1, "lat"].values

        gt3["xc"] = gt3[["x1", "x2", "x3"]].mean(axis=1)  # mean lon of the element
        gt3["yc"] = gt3[["y1", "y2", "y3"]].mean(axis=1)

        # take care of quads
        gt3["x4"] = np.nan
        gt3["y4"] = np.nan

        ide = gt3.loc[gt3.type == 4].index.values  # get the index of finite values

        gt3.loc[ide, "x4"] = grd.loc[gt3.loc[ide, "gd"].values - 1, "lon"].values
        gt3.loc[ide, "y4"] = grd.loc[gt3.loc[ide, "gd"].values - 1, "lat"].values

        gt3.loc[ide, "xc"] = gt3[["x1", "x2", "x3", "x4"]].mean(axis=1)  # mean lon of the element
        gt3.loc[ide, "yc"] = gt3[["y1", "y2", "y3", "y4"]].mean(axis=1)

        ## min kbe
        gt3["kbe4"] = np.nan  # for quads

        gt3["kbe1"] = grd.loc[gt3["ga"] - 1, "kbp00"].values
        gt3["kbe2"] = grd.loc[gt3["gb"] - 1, "kbp00"].values
        gt3["kbe3"] = grd.loc[gt3["gc"] - 1, "kbp00"].values
        gt3.loc[ide, "kbe4"] = grd.loc[gt3.loc[ide, "gd"].values - 1, "kbp00"].values

        gt3["kbe"] = gt3[["kbe1", "kbe2", "kbe3"]].min(axis=1)
        gt3.loc[ide, "kbe"] = gt3[["kbe1", "kbe2", "kbe3", "kbe4"]].min(axis=1)

        self.misc.update({"gt3": gt3.set_index("index")})  # set index back

        # Droping duplicates
        self.misc.update({"melems": elems.loc[elems.global_n.drop_duplicates().index]})  # create the retaining mask
        self.misc.update({"msides": re.loc[re.global_n.drop_duplicates().index]})  # keep only one of the duplicates
        self.misc.update(
            {"mnodes": nodes.loc[nodes.global_n.drop_duplicates().index]}
        )  # keep only one of the duplicates

    def hotstart(self, it=None, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")

        if not "melems" in self.misc:
            self.global2local(**kwargs)

        hfiles = glob.glob(path + "outputs/hotstart_*_{}.nc".format(it))
        hfiles.sort()

        # store them in a list
        out = []
        for i in range(len(hfiles)):
            out.append(xr.open_dataset(hfiles[i]))

        # Create dataset
        side = []
        node = []
        el = []
        one = []

        for key in out[0].variables:
            if "nResident_side" in out[0][key].dims:
                r = self.combine_(key, out, self.misc["msides"], "nResident_side")
                side.append(r)
            elif "nResident_node" in out[0][key].dims:
                r = self.combine_(key, out, self.misc["mnodes"], "nResident_node")
                node.append(r)
            elif "nResident_elem" in out[0][key].dims:
                r = self.combine_(key, out, self.misc["melems"], "nResident_elem")
                el.append(r)
            elif len(out[0][key].dims) == 1:
                one.append(out[0][key])

        side = xr.merge(side).rename({"nResident_side": "side"})
        el = xr.merge(el).rename({"nResident_elem": "elem"})
        node = xr.merge(node).rename({"nResident_node": "node"})
        one = xr.merge(one).rename({"one": "one_new", "it": "iths"})

        # merge
        xdat = xr.merge([side, el, node, one])

        hfile = "hotstart_it={}.nc".format(xdat.iths.values[0])
        logger.info("saving hotstart file\n")

        xdat.to_netcdf(path + "outputs/{}".format(hfile))

    ## Any variable
    def combine(self, out, g2l, name):
        keys = g2l.index.get_level_values(0).unique()
        r = []
        for i in range(len(out)):
            v = out[i].to_pandas()
            v.index += 1
            mask = g2l.loc[keys[i], "local"]
            vm = v.loc[mask]
            vm.index = g2l.loc[keys[i], "global_n"].values
            r.append(vm)
        r = pd.concat(r).sort_index()
        r.index -= 1
        r.index.name = name
        return r

    def combine_(self, var, out, g2l, name):
        if len(out[0][var].shape) == 3:
            dd = []
            for i in range(out[0][var].shape[2]):
                o = []
                for j in range(len(out)):
                    o.append(out[j][var].loc[{out[j][var].dims[2]: i}])
                r = self.combine(o, g2l, name)
                dd.append(xr.DataArray(r.values, dims=out[0][var].dims[:2], name=var))

            tr = xr.concat(dd, dim=out[0][var].dims[2])
            a, b, c = out[0][var].dims
            tr = tr.transpose(a, b, c)
            return tr

        elif len(out[0][var].shape) < 3:
            o = []
            for j in range(len(out)):
                o.append(out[j][var])
            r = self.combine(o, g2l, name)
            return xr.DataArray(r.values, dims=list(out[0][var].dims), name=var)

    def xcombine(self, tfs):

        # Create dataset
        side = []
        node = []
        el = []
        single = []

        for key in tfs[0].variables:
            if "nSCHISM_hgrid_face" in tfs[0][key].dims:
                r = self.combine_(key, tfs, self.misc["melems"], "nSCHISM_hgrid_face")
                side.append(r)
            elif "nSCHISM_hgrid_node" in tfs[0][key].dims:
                r = self.combine_(key, tfs, self.misc["mnodes"], "nSCHISM_hgrid_node")
                node.append(r)
            elif "nSCHISM_hgrid_edge" in tfs[0][key].dims:
                r = self.combine_(key, tfs, self.misc["msides"], "nSCHISM_hgrid_edge")
                node.append(r)
            elif len(tfs[0][key].dims) == 1:
                single.append(tfs[0][key])

        side = xr.merge(side)
        el = xr.merge(el)
        node = xr.merge(node)
        single = xr.merge(single)

        # merge
        return xr.merge([side, el, node, single])

    def tcombine(self, hfiles, sdate, times):

        xall = []
        for k in range(times.shape[0]):
            tfs = []
            for i in range(len(hfiles)):
                tfs.append(xr.open_dataset(hfiles[i]).isel(time=k))

            xall.append(self.xcombine(tfs))

        xdat = xr.concat(xall, dim="time")
        xdat = xdat.assign_coords(time=times)

        return xdat

    # https://stackoverflow.com/questions/41164630/pythonic-way-of-removing-reversed-duplicates-in-list
    @staticmethod
    def remove_reversed_duplicates(iterable):
        # Create a set for already seen elements
        seen = set()
        for item in iterable:
            # Lists are mutable so we need tuples for the set-operations.
            tup = tuple(item)
            if tup not in seen:
                # If the tuple is not in the set append it in REVERSED order.
                seen.add(tup[::-1])
                # If you also want to remove normal duplicates uncomment the next line
                # seen.add(tup)
                yield item

    def read_vgrid(self, **kwargs):

        logger.info("Read vgrid.in\n")

        path = get_value(self, kwargs, "rpath", "./schism/")

        vgrid = pd.read_csv(path + "vgrid.in", header=None, index_col=False, engine="python", delimiter="!")

        try:
            vgrid = vgrid.drop(1, axis=1)
        except:
            pass

        self.misc.update({"ivcor": int(vgrid.iloc[0][0])})

        [Nz, kz, hs] = vgrid.iloc[1].str.split()[0]

        self.misc.update({"Nz": int(Nz)})
        self.misc.update({"kz": int(kz)})
        self.misc.update({"hs": float(hs)})

        zlevels = vgrid.iloc[3 : 3 + self.misc["kz"], 0].str.split(n=2, expand=True)
        zlevels.columns = ["level_index", "z-coordinates"]
        zlevels.set_index("level_index", inplace=True)

        self.misc.update({"zlevels": zlevels})

        constants_index = 3 + self.misc["kz"] + 1

        [h_c, theta_b, theta_f] = vgrid.iloc[constants_index].str.split()[0]
        self.misc.update({"h_c": float(h_c)})
        self.misc.update({"theta_b": float(theta_b)})
        self.misc.update({"theta_f": float(theta_f)})

        sl_index0 = constants_index + 1
        sl_index1 = sl_index0 + self.misc["Nz"] - self.misc["kz"] + 1

        slevels = vgrid.iloc[sl_index0:sl_index1, 0].str.split(n=2, expand=True)
        slevels.columns = ["level_index", "s-coordinates"]
        slevels.set_index("level_index", inplace=True)

        self.misc.update({"slevels": slevels})

    def results(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")

        if len(self.misc) == 0:
            logger.info("retrieving index references ... \n")
            self.global2local(**kwargs)
            logger.info("... done \n")

        # Create grid xarray Dataset
        grd = self.misc["grd"]
        gt3 = self.misc["gt3"]

        # node based variables
        grd.kbp00 = grd.kbp00.astype(int)
        xnodes = grd.to_xarray().rename(
            {
                "lon": "SCHISM_hgrid_node_x",
                "lat": "SCHISM_hgrid_node_y",
                "kbp00": "node_bottom_index",
                "index": "nSCHISM_hgrid_node",
            }
        )

        xnodes = xnodes.drop_vars("nSCHISM_hgrid_node")

        # element based variables
        gt34 = gt3.loc[:, ["ga", "gb", "gc", "gd"]].values  # SCHISM_hgrid_face_nodes
        xelems = xr.Dataset(
            {
                "SCHISM_hgrid_face_nodes": (
                    ["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
                    gt34 - 1,
                ),  # -> start index = 0
                "SCHISM_hgrid_face_x": (["nSCHISM_hgrid_face"], gt3.loc[:, "xc"].values),
                "SCHISM_hgrid_face_y": (["nSCHISM_hgrid_face"], gt3.loc[:, "yc"].values),
                "ele_bottom_index": (["nSCHISM_hgrid_face"], gt3.kbe.values),
            }
        )

        logger.info("done with node based variables \n")

        # edge based variables
        sides = []
        for [etype, ga, gb, gc, gd] in gt3.loc[:, ["type", "ga", "gb", "gc", "gd"]].values:
            if etype == 3:
                sides.append([gb, gc])
                sides.append([gc, ga])
                sides.append([ga, gb])
            elif etype == 4:
                sides.append([gb, gc])
                sides.append([gc, gd])
                sides.append([gd, ga])
                sides.append([ga, gb])

        # removing duplicates
        sides = list(self.remove_reversed_duplicates(sides))

        ed = pd.DataFrame(sides, columns=["node1", "node2"])

        # mean x, y
        ed["x1"] = grd.loc[ed["node1"].values - 1, "lon"].values  # lon of the index, -1 for python convention
        ed["y1"] = grd.loc[ed["node1"].values - 1, "lat"].values  # lat of the index
        ed["x2"] = grd.loc[ed["node2"].values - 1, "lon"].values
        ed["y2"] = grd.loc[ed["node2"].values - 1, "lat"].values

        ed["xc"] = ed[["x1", "x2"]].mean(axis=1)  # mean of the edge index
        ed["yc"] = ed[["y1", "y2"]].mean(axis=1)

        ## min bottom index
        ed["kbs1"] = grd.loc[ed["node1"] - 1, "kbp00"].values
        ed["kbs2"] = grd.loc[ed["node2"] - 1, "kbp00"].values

        ed["kbs"] = ed[["kbs1", "kbs2"]].min(axis=1)

        xsides = xr.Dataset(
            {
                "SCHISM_hgrid_edge_nodes": (
                    ["nSCHISM_hgrid_edge", "two"],
                    [[x - 1, y - 1] for [x, y] in sides],
                ),  # index from 0
                "SCHISM_hgrid_edge_x": (["nSCHISM_hgrid_edge"], ed["xc"].values),
                "SCHISM_hgrid_edge_y": (["nSCHISM_hgrid_edge"], ed["yc"].values),
                "edge_bottom_index": (["nSCHISM_hgrid_edge"], ed.kbs.values),
            }
        )

        logger.info("done with side based variables \n")

        # General properties

        header2 = self.misc["header"].apply(pd.to_numeric)
        nlist = [
            "start_year",
            "start_month",
            "start_day",
            "start_hour",
            "utc_start",
            "dtout",
            "nspool",
            "nvrt",
            "kz",
            "ics",
        ]
        header2[nlist] = header2[nlist].astype(int)
        sigmas = [x for x in header2.columns if "sigma" in x]
        sigms = header2.loc[:, sigmas].values.flatten()  # get sigmas
        iwet_dry = 0  # defined by the user
        ihgrid_id = -2147483647  # defined by user - 0,dummy_dim,ihgrid_id
        one = xr.Dataset({"dry_value_flag": (("one"), [iwet_dry]), "SCHISM_hgrid": (("one"), [ihgrid_id])})

        # compute cs
        klev = np.arange(header2.kz.values[0], header2.nvrt.values[0] + 1)
        k = klev - header2.kz.values

        cs = np.zeros(k)

        cs = (1 - header2.theta_b.values) * np.sinh(header2.theta_f.values * sigms[k]) / np.sinh(
            header2.theta_f.values
        ) + header2.theta_b.values * (
            np.tanh(header2.theta_f.values * (sigms[k] + 0.5)) - np.tanh(header2.theta_f.values * 0.5)
        ) / 2 / np.tanh(
            header2.theta_f.values * 0.5
        )

        Cs = xr.Dataset({"Cs": (("sigma"), cs)}, coords={"sigma": sigms})

        header_list = ["ics", "h0", "h_c", "theta_b", "theta_f", "h_s"]
        gen = header2[header_list].to_xarray()
        gen = gen.rename({"index": "one"})

        # merge
        gen = xr.merge([gen, Cs, one])

        gen = gen.rename(
            {
                "ics": "coordinate_system_flag",
                "h0": "minimum_depth",
                "h_c": "sigma_h_c",
                "theta_b": "sigma_theta_b",
                "theta_f": "sigma_theta_f",
                "h_s": "sigma_maxdepth",
            }
        )

        gen = gen.drop_vars("one")

        # set timestamp
        date = header2.loc[:, ["start_year", "start_month", "start_day", "start_hour", "utc_start"]]
        date = date.astype(int)
        date.columns = ["year", "month", "day", "hour", "utc"]  # rename the columns
        # set the start timestamp
        sdate = pd.Timestamp(
            year=date.year.values[0],
            month=date.month.values[0],
            day=date.day.values[0],
            hour=date.hour.values[0],
            tz=date.utc.values[0],
        )

        logger.info("done with generic variables \n")

        # Read Netcdf output files
        hfiles = glob.glob(path + "outputs/schout_*_*.nc")

        irange_ = [int(x.split("_")[-1].split(".")[0]) for x in hfiles]
        irange_ = np.unique(irange_)

        irange = get_value(self, kwargs, "rlist", irange_)

        self.read_vgrid()  # read grid attributes

        logger.info("Write combined NetCDF files \n")
        #        total_xdat = []
        for val in tqdm(irange):
            hfiles = glob.glob(path + "outputs/schout_*_{}.nc".format(val))
            hfiles.sort()

            times = xr.open_dataset(hfiles[0]).time
            times = pd.to_datetime(times.values, unit="s", origin=sdate.tz_convert(None))

            if times.size == 0:
                continue

            idat = self.tcombine(hfiles, sdate, times)

            # MERGE

            xc = xr.merge([idat, gen, xnodes, xelems, xsides])

            # Choose attrs
            if header2.ics.values == 1:
                lat_coord_standard_name = "projection_y_coordinate"
                lon_coord_standard_name = "projection_x_coordinate"
                x_units = "m"
                y_units = "m"
                lat_str_len = 23
                lon_str_len = 23
            else:
                lat_coord_standard_name = "latitude"
                lon_coord_standard_name = "longitude"
                x_units = "degrees_east"
                y_units = "degrees_north"
                lat_str_len = 8
                lon_str_len = 9

            # set Attrs
            xc.SCHISM_hgrid_node_x.attrs = {
                "long_name": "node x-coordinate",
                "standard_name": lon_coord_standard_name,
                "units": x_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.SCHISM_hgrid_node_y.attrs = {
                "long_name": "node y-coordinate",
                "standard_name": lat_coord_standard_name,
                "units": y_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.depth.attrs = {
                "long_name": "Bathymetry",
                "units": "meters",
                "positive": "down",
                "mesh": "SCHISM_hgrid",
                "location": "node",
            }

            xc.sigma_h_c.attrs = {
                "long_name": "ocean_s_coordinate h_c constant",
                "units": "meters",
                "positive": "down",
            }

            xc.sigma_theta_b.attrs = {"long_name": "ocean_s_coordinate theta_b constant"}

            xc.sigma_theta_f.attrs = {"long_name": "ocean_s_coordinate theta_f constant"}

            xc.sigma_maxdepth.attrs = {
                "long_name": "ocean_s_coordinate maximum depth cutoff (mixed s over z boundary)",
                "units": "meters",
                "positive": "down",
            }

            xc.Cs.attrs = {"long_name": "Function C(s) at whole levels", "positive": "up"}

            xc.dry_value_flag.attrs = {"values": "0: use last-wet value; 1: use junk"}

            xc.SCHISM_hgrid_face_nodes.attrs = {
                "long_name": "Horizontal Element Table",
                "cf_role": "face_node_connectivity",
                "start_index": 0,
            }

            xc.SCHISM_hgrid_edge_nodes.attrs = {
                "long_name": "Map every edge to the two nodes that it connects",
                "cf_role": "edge_node_connectivity",
                "start_index": 0,
            }

            xc.SCHISM_hgrid_edge_x.attrs = {
                "long_name": "x_coordinate of 2D mesh edge",
                "standard_name": lon_coord_standard_name,
                "units": x_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.SCHISM_hgrid_edge_y.attrs = {
                "long_name": "y_coordinate of 2D mesh edge",
                "standard_name": lat_coord_standard_name,
                "units": y_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.SCHISM_hgrid_face_x.attrs = {
                "long_name": "x_coordinate of 2D mesh face",
                "standard_name": lon_coord_standard_name,
                "units": x_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.SCHISM_hgrid_face_y.attrs = {
                "long_name": "y_coordinate of 2D mesh face",
                "standard_name": lat_coord_standard_name,
                "units": y_units,
                "mesh": "SCHISM_hgrid",
            }

            xc.SCHISM_hgrid.attrs = {
                "long_name": "Topology data of 2d unstructured mesh",
                "topology_dimension": 2,
                "cf_role": "mesh_topology",
                "node_coordinates": "SCHISM_hgrid_node_x SCHISM_hgrid_node_y",
                "face_node_connectivity": "SCHISM_hgrid_face_nodes",
                "edge_coordinates": "SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y",
                "face_coordinates": "SCHISM_hgrid_face_x SCHISM_hgrid_face_y",
                "edge_node_connectivity": "SCHISM_hgrid_edge_nodes",
            }

            xc.node_bottom_index.attrs = {
                "long_name": "bottom level index at each node",
                "units": "non-dimensional",
                "mesh": "SCHISM_hgrid",
                "location": "node",
                "start_index": 0,
            }

            xc.ele_bottom_index.attrs = {
                "long_name": "bottom level index at each element",
                "units": "non-dimensional",
                "mesh": "SCHISM_hgrid",
                "location": "elem",
                "start_index": 0,
            }

            xc.edge_bottom_index.attrs = {
                "long_name": "bottom level index at each edge",
                "units": "non-dimensional",
                "mesh": "SCHISM_hgrid",
                "location": "edge",
                "start_index": 0,
            }

            base_date = " ".join([str(x) for x in date.T.values.flatten()])
            xc.time.attrs = {"long_name": "Time", "base_date": base_date, "standard_name": "time"}

            xc.sigma.attrs = {
                "long_name": "S coordinates at whole levels",
                "units": "1",
                "standard_name": "ocean_s_coordinate",
                "positive": "up",
                "h_s": self.misc["hs"],
                "h_c": self.misc["h_c"],
                "theta_b": self.misc["theta_b"],
                "theta_f": self.misc["theta_f"],
                "formula_terms": "s: sigma eta: elev depth: depth a: sigma_theta_f b: sigma_theta_b depth_c: sigma_h_c",
            }

            # Dataset Attrs

            xc.attrs = {
                "Conventions": "CF-1.0, UGRID-1.0",
                "title": "SCHISM Model output",
                "source": "SCHISM model output version v10",
                "references": "http://ccrm.vims.edu/schismweb/",
                "history": "created by pyposeidon",
                "comment": "SCHISM Model output",
                "type": "SCHISM Model output",
                "VisIT_plugin": "https://schism.water.ca.gov/library/-/document_library/view/3476283",
            }

            xc.to_netcdf(path + "outputs/schout_{}.nc".format(val))

        logger.info("done with output netCDF files \n")

    def set_obs(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")
        nspool_sta = get_value(self, kwargs, "nspool_sta", 1)
        tg_database = get_value(self, kwargs, "obs", DATA_PATH + "critech.csv")
        coastal_monitoring = get_value(self, kwargs, "coastal_monitoring", False)
        flags = get_value(self, kwargs, "station_flags", [1] + [0] * 8)

        station_flag = pd.DataFrame(
            {
                "elev": flags[0],
                "air_pressure": flags[1],
                "windx": flags[2],
                "windy": flags[3],
                "T": flags[4],
                "S": flags[5],
                "u": flags[6],
                "v": flags[7],
                "w": flags[8],
            },
            index=[0],
        )

        ## FOR TIDE GAUGE MONITORING
        z = self.__dict__.copy()
        tg = obs(**z)

        logger.info("get in-situ measurements locations \n")

        gpoints = np.array(
            list(zip(self.grid.Dataset.SCHISM_hgrid_node_x.values, self.grid.Dataset.SCHISM_hgrid_node_y.values))
        )

        stations = []
        grid_index = []
        for l in range(tg.locations.shape[0]):
            plat, plon = tg.locations.loc[l, ["latitude", "longitude"]]
            cp = closest_node([plon, plat], gpoints)
            grid_index.append(
                list(
                    zip(self.grid.Dataset.SCHISM_hgrid_node_x.values, self.grid.Dataset.SCHISM_hgrid_node_y.values)
                ).index(tuple(cp))
            )
            stations.append(cp)

        # to df
        stations = pd.DataFrame(stations, columns=["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y"])
        stations["z"] = 0
        stations.index += 1
        stations["gindex"] = grid_index

        if coastal_monitoring:
            ## FOR COASTAL MONITORING
            logger.info("set land boundaries as observation points \n")

            # get land boundaries
            coasts = [key for key in self.grid.Dataset.variables if "land_boundary" in key]
            # get index
            bnodes = []
            for i in range(len(coasts)):
                bnodes = bnodes + list(self.grid.Dataset[coasts[i]].dropna("index").astype(int).values)
            bnodes = np.unique(bnodes)  # keep unique values
            # get x,y
            xx = self.grid.Dataset.SCHISM_hgrid_node_x.to_dataframe().loc[bnodes]
            yy = self.grid.Dataset.SCHISM_hgrid_node_y.to_dataframe().loc[bnodes]
            # put them in dataframe
            coastal_stations = pd.concat([xx, yy], axis=1, sort=False)
            coastal_stations["z"] = 0
            coastal_stations.reset_index(inplace=True, drop=True)
            coastal_stations.index += 1
            coastal_stations["gindex"] = bnodes
            coastal_stations.head()
            # append
            stations = pd.concat([stations, coastal_stations])
            stations.reset_index(inplace=True, drop=True)
            stations.index += 1

        stations["gindex"] = stations["gindex"].astype(int)

        # modify config paramater
        self.params["SCHOUT"]["iout_sta"] = 1
        self.params["SCHOUT"]["nspool_sta"] = nspool_sta
        self.params.write(path + "param.nml", force=True)

        self.stations = stations

        logger.info("write out stations.in file \n")

        # output to file
        with open(path + "station.in", "w") as f:
            station_flag.to_csv(f, header=None, index=False, sep=" ")
            f.write("{}\n".format(stations.shape[0]))
            stations.loc[:, ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "z"]].to_csv(f, header=None, sep=" ")

    def get_station_data(self, **kwargs):

        path = get_value(self, kwargs, "rpath", "./schism/")

        # locate the station files
        sfiles = glob.glob(path + "outputs/staout_*")
        sfiles.sort()

        try:
            # get the station flags
            flags = pd.read_csv(path + "station.in", header=None, nrows=1, delim_whitespace=True).T
            flags.columns = ["flag"]
            flags["variable"] = ["elev", "air_pressure", "windx", "windy", "T", "S", "u", "v", "w"]

            vals = flags[flags.values == 1]  # get the active ones
        except OSError as e:
            if e.errno == errno.EEXIST:
                logger.error("No station.in file present")
            return

        dstamp = kwargs.get("dstamp", self.start_date)

        dfs = []
        for idx in vals.index:
            obs = np.loadtxt(sfiles[idx])
            df = pd.DataFrame(obs)
            df = df.set_index(0)
            df.index.name = "time"
            df.columns.name = vals.loc[idx, "variable"]
            df.index = pd.to_datetime(dstamp) + pd.to_timedelta(df.index, unit="S")
            pindex = pd.MultiIndex.from_product([df.T.columns, df.T.index])

            r = pd.DataFrame(df.values.flatten(), index=pindex, columns=[vals.loc[idx, "variable"]])
            r.index.names = ["time", "index"]

            dfs.append(r.to_xarray())

        self.time_series = xr.combine_by_coords(dfs)

    def get_data(self, **kwargs):

        dic = self.__dict__

        dic.update(kwargs)

        self.data = data(**dic)
