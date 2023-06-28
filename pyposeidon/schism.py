"""
Schism model of pyposeidon. It controls the creation, execution & output  of a complete simulation based on SCHISM.

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
import shlex
import sys
import json
from collections import OrderedDict
import pandas as pd
import geopandas as gp
import glob
from shutil import copyfile
import xarray as xr
import geopandas as gp
import shapely
import f90nml
import errno
import dask
from searvey import ioc
from tqdm.auto import tqdm

# local modules
import pyposeidon
import pyposeidon.mesh as pmesh
import pyposeidon.meteo as pmeteo
import pyposeidon.dem as pdem
from pyposeidon.paths import DATA_PATH
from pyposeidon.utils.get_value import get_value
from pyposeidon.utils.converter import myconverter
from pyposeidon.utils.cpoint import closest_node
from pyposeidon.utils.unml import unml
from pyposeidon.utils import data
from pyposeidon.utils.norm import normalize_column_names
from pyposeidon.utils.obs import get_obs_data


import logging

from . import tools

logger = logging.getLogger(__name__)


SCHISM_NAME = "schism"


class Schism:
    def __init__(self, **kwargs):
        """
        Create a Schism solver

        !!! danger ""
            Due to a limitation of the Library rendering the docstrings, all arguments are marked
            as `required`, nevertheless they are all `Optional`.

        Args:
            rfolder str: The path to a directory containing the results of a Model solved by Schism.
            geometry Union[dict, str, GeoDataFrame]: A `GeoDataFrame` or the path to a shapefile or
                a dict defining the lat/lon window.
            load_mesh bool: Flag indicating whether to load the mesh or not. Defaults to `False`.
            load_meteo bool: Flag indicating whether to load the meteo data or not. Defauls to
                `False`.
            coastlines Union[str, GeoDataFrame]: A `GeoDataFrame` or the path to a shapefile which
                describes the coastlines.
            tag str: The model's "tag". Defaults to `"schism"`.
            tide str: Flag indicating whether to load "tide". Defaults to `False`.
            atm bool: The solver's atm. Defaults to `True`.
            monitor bool: The solver's monitor. Defaults to `False`.
            epath str: The path to the schism executable. If the `SCHISM` env variable has been
                set, then it overrides the value passed as the parameter.
            start_date str: The date from which the analysis should start. It should be a string parseable
                by `pd.to_datetime()`.
            end_date str: The date at which the analysis should end. It should be a string parseable by
                `pd.to_datetime()`.
            time_frame str: The duration of the analysis. It should be a string parseable by
                `pd.to_datetime()`.
            rdate str: Reference date of the run.
            rpath str: Path for output of the model. Defaults to `./schism/`.
            m_index int: Define the index of the meteo Dataset. Defaults to `1`.
            filename str: Path to output the meteo Dataset. Defaults to `sflux/`.
            dstamp str: Reference date for station data. Defaults to rdate.
            parameters dict: Overwrite default Schism's parameter values.
            meteo_source str: Path or url to meteo data.
            dem_source str: Path or url to bathymetric data.
            update list[str]: Control the update of the model e.g `['dem']`-> updates only bathymetry.
                Defaults to `["all"]`.
            meteo_split_by str: Split the meteo Dataset to multiple files by e.g. `"day"`.
                Defaults to `None`.
            manning_file str: Path to manning file.
            manning float: Set constant value in the manning.gr3 file. Defaults to `0.12`.
            windrot_file str: Path to windrot file.
            windrot float: Set constant value in the windrot_geo2proj.gr3 file. Defaults to `0.00001`.
            station_flags list[int]: Define the flag for station output. Defaults to `[1,0,0,0,0,0,0,0,0]`.
            coastal_monitoring bool: Flag for setting all land/island boundary nodes as stations.
                Defaults to `False`.
            obs str: Path to csv file for station locations. Defaults to `searvey`.
            nspool_sta int: Related to station nodes setup. Defaults to `1`.
        """

        rfolder = kwargs.get("rfolder", None)
        if rfolder:
            self.read_folder(**kwargs)

        self.geometry = kwargs.get("geometry", None)

        if self.geometry:
            if isinstance(self.geometry, dict):
                self.lon_min = self.geometry["lon_min"]
                self.lon_max = self.geometry["lon_max"]
                self.lat_min = self.geometry["lat_min"]
                self.lat_max = self.geometry["lat_max"]
            elif self.geometry == "global":
                logger.warning("geometry is 'global'")
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
            else:
                logger.warning("no geometry given")

        # coastlines
        coastlines = kwargs.get("coastlines", None)

        if coastlines is not None:
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

        if not hasattr(self, "rdate"):
            self.rdate = get_value(self, kwargs, "rdate", self.start_date)

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

        self.solver_name = SCHISM_NAME

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

            config_file = os.path.join(DATA_PATH, "param.nml")

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
            self.parameters.update(dic)

        # test rnday
        if float(params["CORE"]["rnday"]) * 24 * 3600 > (self.end_date - self.start_date).total_seconds():
            # ---------------------------------------------------------------------
            logger.warning("rnday larger than simulation range\n")
            logger.warning(
                "rnday={} while simulation time is {}\n".format(
                    params["core"]["rnday"],
                    (self.end_date - self.start_date).total_seconds() / (3600 * 24.0),
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
            self.params.write(os.path.join(path, "param.nml"), force=True)

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
                self.meteo = pmeteo.Meteo(**z)
            else:
                logger.info("skipping meteo ..\n")
        else:
            self.meteo = pmeteo.Meteo(**z)

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

        zero = dask.array.zeros(ar[p].data.shape)

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

        sout.attrs = {
            "description": "Schism meteo data",
            "history": "pyposeidon",
            "source": "netCDF4 python module",
        }

        sout.time.attrs = {
            "long_name": "Time",
            "standard_name": "time",
            "base_date": bdate,
            "units": udate,
        }

        sout.lat.attrs = {
            "units": "degrees_north",
            "long_name": "Latitude",
            "standard_name": "latitude",
        }

        sout.lon.attrs = {
            "units": "degrees_east",
            "long_name": "Longitude",
            "standard_name": "longitude",
        }

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

        # Create sflux directory if necessary
        sflux_path = os.path.join(path, "sflux")
        os.makedirs(sflux_path, exist_ok=True)

        # Save meteo netcdf to disk
        m_index = kwargs.get("m_index", 1)
        filename = kwargs.get("filename", f"sflux_air_{m_index}.0001.nc")
        netcdf_path = os.path.join(sflux_path, filename)
        sout.to_netcdf(netcdf_path)

    # ============================================================================================
    # DEM
    # ============================================================================================

    def bath(self, **kwargs):
        #       z = self.__dict__.copy()

        if self.mesh.Dataset is not None:
            kwargs["grid_x"] = self.mesh.Dataset.SCHISM_hgrid_node_x.values
            kwargs["grid_y"] = self.mesh.Dataset.SCHISM_hgrid_node_y.values

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
                self.dem.Dataset = pdem.dem_on_mesh(self.dem.Dataset, **kwargs)

                if kwargs.get("adjust_dem", True):
                    coastline = kwargs.get("coastlines", None)
                    if coastline is None:
                        logger.warning("coastlines not present, aborting adjusting dem\n")
                    elif "adjusted" in self.dem.Dataset.data_vars.keys():
                        logger.info("Dem already adjusted\n")
                    elif "fval" in self.dem.Dataset.data_vars.keys():
                        logger.info("Dem already adjusted\n")
                    else:
                        self.dem.adjust(coastline, **kwargs)
                try:
                    try:
                        bat = -self.dem.Dataset.fval.values.astype(float)  # minus for the hydro run
                        if np.isnan(bat).sum() != 0:
                            raise Exception("Bathymetry contains NaNs")
                            logger.warning("Bathymetric values fval contain NaNs, using ival values ..\n")

                    except:
                        bat = -self.dem.Dataset.ival.values.astype(float)  # minus for the hydro run

                    self.mesh.Dataset.depth.loc[: bat.size] = bat

                    logger.info("updating bathymetry ..\n")

                except AttributeError as e:
                    logger.info("Keeping bathymetry in hgrid.gr3 due to {}\n".format(e))

            else:
                logger.info("dem from mesh file\n")

    # ============================================================================================
    # EXECUTION
    # ============================================================================================
    def create(self, **kwargs):
        if not kwargs:
            kwargs = self.__dict__.copy()

        # Set background dem as scale for mesh generation
        dpath = get_value(self, kwargs, "dem_source", None)

        if dpath:
            self.dem = pdem.Dem(**kwargs)
            kwargs.update({"dem_source": self.dem.Dataset})
        else:
            logger.info("no dem available\n")

        # Mesh
        self.mesh = pmesh.set(type="tri2d", **kwargs)

        # set lat/lon from file
        if self.mesh.Dataset is not None:
            kwargs.update({"lon_min": self.mesh.Dataset.SCHISM_hgrid_node_x.values.min()})
            kwargs.update({"lon_max": self.mesh.Dataset.SCHISM_hgrid_node_x.values.max()})
            kwargs.update({"lat_min": self.mesh.Dataset.SCHISM_hgrid_node_y.values.min()})
            kwargs.update({"lat_max": self.mesh.Dataset.SCHISM_hgrid_node_y.values.max()})

            self.lon_min = self.mesh.Dataset.SCHISM_hgrid_node_x.values.min()
            self.lon_max = self.mesh.Dataset.SCHISM_hgrid_node_x.values.max()
            self.lat_min = self.mesh.Dataset.SCHISM_hgrid_node_y.values.min()
            self.lat_max = self.mesh.Dataset.SCHISM_hgrid_node_y.values.max()

        # get bathymetry
        if self.mesh.Dataset is not None:
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
        sflux_path = os.path.join(path, "sflux")
        os.makedirs(sflux_path, exist_ok=True)

        with open(os.path.join(sflux_path, "sflux_inputs.txt"), "w") as f:
            f.write("&sflux_inputs\n")
            f.write("/ \n\n")

        # save params.in

        self.params.write(os.path.join(path, "param.nml"), force=True)

        # Mesh related files
        if self.mesh.Dataset is not None:
            # save bctides.in
            bs = self.mesh.Dataset[["node", "id", "type"]].to_dataframe()
            # open boundaries
            number_of_open_boundaries = bs.loc[bs.type == "open"].id
            if not number_of_open_boundaries.empty:
                number_of_open_boundaries = number_of_open_boundaries.max()
            else:
                number_of_open_boundaries = 0
            number_of_open_boundaries_nodes = bs.loc[bs.type == "open"].shape[0]

            with open(os.path.join(path, "bctides.in"), "w") as f:
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
            with open(os.path.join(path, "vgrid.in"), "w") as f:
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

            # save hgrid.gr3
            try:
                try:
                    bat = -self.dem.Dataset.fval.values.astype(float)  # minus for the hydro run
                    if np.isnan(bat).sum() != 0:
                        raise Exception("Bathymetry contains NaNs")
                        logger.warning("Bathymetric values fval contain NaNs, using ival values ..\n")

                except:
                    bat = -self.dem.Dataset.ival.values.astype(float)  # minus for the hydro run

                self.mesh.Dataset.depth.loc[: bat.size] = bat

                self.mesh.to_file(filename=os.path.join(path, "hgrid.gr3"))
                copyfile(os.path.join(path, "hgrid.gr3"), os.path.join(path, "hgrid.ll"))

                logger.info("updating bathymetry ..\n")

            except AttributeError as e:
                logger.info("Keeping bathymetry from hgrid.gr3 ..\n")

                copyfile(self.mesh_file, os.path.join(path, "hgrid.gr3"))  # copy original grid file
                copyfile(os.path.join(path, "hgrid.gr3"), os.path.join(path, "hgrid.ll"))

            # manning file
            manfile = os.path.join(path, "manning.gr3")

            if hasattr(self, "manning_file"):
                copyfile(self.manning_file, manfile)  # copy original manning file
                if self.manning_file == manfile:
                    logger.info("Keeping manning file ..\n")

            manning = get_value(self, kwargs, "manning", 0.12)
            nn = self.mesh.Dataset.nSCHISM_hgrid_node.size
            n3e = self.mesh.Dataset.nSCHISM_hgrid_face.size

            with open(manfile, "w") as f:
                f.write("\t 0 \n")
                f.write("\t {} {}\n".format(n3e, nn))

            df = self.mesh.Dataset[["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]].to_dataframe()

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

            windfile = os.path.join(path, "windrot_geo2proj.gr3")

            if hasattr(self, "windrot_file"):
                copyfile(self.windrot_file, windfile)  # copy original grid file
                if self.windrot_file != windfile:
                    logger.info("Keeping windrot_geo2proj file ..\n")

            windrot = get_value(self, kwargs, "windrot", 0.00001)

            with open(windfile, "w") as f:
                f.write("\t 0 \n")
                f.write("\t {} {}\n".format(n3e, nn))

            df = self.mesh.Dataset[["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]].to_dataframe()

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
                    mpaths = ["sflux_air_{}.{:04d}.nc".format(m_index, t + 1) for t in np.arange(len(times))]
                    for das, mpath in list(zip(datasets, mpaths)):
                        self.to_force(
                            das, vars=["msl", "u10", "v10"], rpath=path, filename=mpath, date=self.rdate, **kwargs
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
            logger.warning("schism executable path (epath) not given -> using default \n")
            # ------------------------------------------------------------------------------
            bin_path = "schism"

        scribes = get_value(self, kwargs, "scribes", 0)

        tools.create_schism_mpirun_script(
            target_dir=calc_dir, cmd=bin_path, script_name="launchSchism.sh", scribes=scribes
        )

        # ---------------------------------------------------------------------
        logger.info("output done\n")
        # ---------------------------------------------------------------------

    def run(self, **kwargs):
        calc_dir = get_value(self, kwargs, "rpath", "./schism/")

        # ---------------------------------------------------------------------
        logger.info("executing model\n")
        # ---------------------------------------------------------------------

        if not tools.is_mpirun_installed():
            logger.warning("mpirun is not installed, ending.. \n")
            return

        cmd = "./launchSchism.sh"

        # note that cwd is the folder where the executable is
        proc = subprocess.run(
            shlex.split(cmd),
            check=False,
            capture_output=True,
            text=True,
            cwd=calc_dir,
        )

        with open(os.path.join(calc_dir, "err.log"), "w") as fd:
            fd.write(proc.stderr)
        with open(os.path.join(calc_dir, "run.log"), "w") as fd:
            fd.write(proc.stdout)

        # store output in class
        self.stderr = proc.stderr
        self.stdout = proc.stdout

        if proc.returncode == 0:
            if ("ABORT" in proc.stderr) or ("ABORT" in proc.stdout):
                # ---------------------------------------------------------------------
                logger.error("schism failed to execute correctly. See logs\n")
                # ---------------------------------------------------------------------
            else:
                # ---------------------------------------------------------------------
                logger.info("finished\n")
                # ---------------------------------------------------------------------
        else:
            # ---------------------------------------------------------------------
            logger.error("schism failed to execute. See logs\n")
            # ---------------------------------------------------------------------
            proc.check_returncode()

    def save(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./schism/")

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
            dic.update({"dem_attrs": dem.Dataset.elevation.attrs})

        meteo = self.__dict__.get("meteo", None)
        if isinstance(meteo, str):
            dic.update({"meteo": meteo})
        elif isinstance(meteo, pmeteo.Meteo):
            try:
                dic.update({"meteo": [meteo.Dataset.attrs]})
            except:
                dic.update({"meteo": [x.attrs for x in meteo.Dataset]})

        coastlines = self.__dict__.get("coastlines", None)
        dic.update({"coastlines": coastlines})

        dic["version"] = pyposeidon.__version__

        for attr, value in dic.items():
            if isinstance(value, datetime.datetime):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.Timedelta):
                dic[attr] = value.isoformat()
            if isinstance(value, pd.DataFrame):
                dic[attr] = value.to_dict()
            if isinstance(value, gp.GeoDataFrame):
                dic[attr] = value.to_json()

        filename = os.path.join(path, f"{self.tag}_model.json")
        json.dump(dic, open(filename, "w"), indent=4, default=myconverter)

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
                logger.warning("schism executable path (epath) not given -> using default \n")
                # ------------------------------------------------------------------------------
                bin_path = "schism"

            tools.create_mpirun_script(
                target_dir=calc_dir,
                cmd=bin_path,
                script_name="launchSchism.sh",
            )

            self.run(**kwargs)

    def read_folder(self, rfolder, **kwargs):
        self.rpath = rfolder
        s = glob.glob(os.path.join(rfolder, "/param.nml"))
        mfiles1 = glob.glob(os.path.join(rfolder, "/sflux/*_1*.nc"))
        mfiles1.sort()
        # check for 2nd meteo
        mfiles2 = glob.glob(os.path.join(rfolder, "/sflux/*_2*.nc"))
        mfiles2.sort()

        mfiles = {"1": mfiles1, "2": mfiles2}

        mfiles = {k: v for k, v in mfiles.items() if v}  # remove empty keys, e.g. no mfiles2

        hfile = os.path.join(rfolder, "/hgrid.gr3")  # Grid
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

        load_mesh = get_value(self, kwargs, "load_mesh", False)

        if load_mesh:
            try:
                self.mesh = pmesh.set(type="tri2d", mesh_file=hfile)
            except:
                logger.warning("loading mesh failed")
                pass
        else:
            logger.warning("no mesh loaded")

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
                self.meteo = pmeteo.Meteo(meteo_source=pm)

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
                            logger.warning("loading meteo failed")
                            break
                    if ma:
                        msource = xr.merge(ma)
                        pm.append(msource)

                if len(pm) == 1:
                    pm = pm[0]

                self.meteo = pmeteo.Meteo(meteo_source=pm)

        else:
            logger.warning("no meteo loaded")

    def global2local(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./schism/")

        # Read the global node index distribution to the cores
        gfiles = glob.glob(os.path.join(path, "outputs/local_to_global_*"))
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
                    f,
                    skiprows=3,
                    header=None,
                    nrows=nels[i],
                    names=["local", "global_n"],
                    delim_whitespace=True,
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
                names=[
                    "start_year",
                    "start_month",
                    "start_day",
                    "start_hour",
                    "utc_start",
                ],
            )
        with open(gfiles[0], "r") as f:
            h1 = pd.read_csv(
                f,
                skiprows=nels[0] + nq[0] + nw[0] + 7,
                header=None,
                nrows=1,
                delim_whitespace=True,
                names=[
                    "nrec",
                    "dtout",
                    "nspool",
                    "nvrt",
                    "kz",
                    "h0",
                    "h_s",
                    "h_c",
                    "theta_b",
                    "theta_f",
                    "ics",
                ],
            )

        ztots = ["ztot_" + str(i) for i in range(1, h1.loc[:, "kz"].values[0] - 1)]
        sigmas = ["sigma_" + str(i) for i in range(h1.loc[:, "nvrt"].values[0] - h1.loc[:, "kz"].values[0] + 1)]

        # read secondary header
        with open(gfiles[0], "r") as f:
            h2 = pd.read_csv(
                f,
                skiprows=nels[0] + nq[0] + nw[0] + 8,
                header=None,
                nrows=2,
                delim_whitespace=True,
            )

        h2 = h2.T
        h2.columns = ztots + sigmas

        # combine headers
        self.misc.update({"header": pd.concat([h0, h1, h2], axis=1)})

        # read lat/lon from all files
        gframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                gframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + nq[i] + nw[i] + 11,
                    header=None,
                    nrows=nq[i],
                    delim_whitespace=True,
                    names=["lon", "lat", "depth", "kbp00"],
                )

        mesh = pd.concat(gframes, keys=keys)

        # Droping duplicates
        drops = nodes.reset_index()[nodes.reset_index().duplicated("global_n")].index.to_list()
        cnodes = nodes.global_n.drop_duplicates()  # drop duplicate global nodes and store the values to an array

        mesh = mesh.reset_index().drop(drops)  # Use the mask from nodes to match mesh
        mesh = mesh.set_index(["level_0", "level_1"])

        mesh.index = mesh.index.droplevel()  # drop multi-index
        mesh = mesh.reset_index(drop=True)  # reset index
        mesh.index = cnodes.values - 1  # reindex based on the global index, -1 for the python convention
        grd = mesh.sort_index()  # sort with the new index (that is the global_n)
        self.misc.update({"grd": grd.reset_index(drop=True)})  # reindex for final version

        # Read tessalation
        eframes = np.empty(len(gfiles), dtype=object)
        for i in range(len(gfiles)):
            with open(gfiles[i], "r") as f:
                eframes[i] = pd.read_csv(
                    f,
                    skiprows=nels[i] + nq[i] + nw[i] + nq[i] + 11,
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

        hfiles = glob.glob(os.path.join(path, f"outputs/hotstart_*_{it}.nc"))
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

        xdat.to_netcdf(os.path.join(path, f"outputs/{hfile}"))

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
                el.append(r)
            elif "nSCHISM_hgrid_node" in tfs[0][key].dims:
                r = self.combine_(key, tfs, self.misc["mnodes"], "nSCHISM_hgrid_node")
                node.append(r)
            elif "nSCHISM_hgrid_edge" in tfs[0][key].dims:
                r = self.combine_(key, tfs, self.misc["msides"], "nSCHISM_hgrid_edge")
                side.append(r)
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
        logger.info("read vgrid.in\n")

        path = get_value(self, kwargs, "rpath", "./schism/")

        vgrid = pd.read_csv(
            os.path.join(path, "vgrid.in"),
            header=None,
            index_col=False,
            engine="python",
            delimiter="!",
        )

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

        logger.info("get combined 2D netcdf files \n")
        hfiles = glob.glob(os.path.join(path, "outputs/out2d_*.nc"))
        hfiles.sort()

        if hfiles:
            x2d = xr.open_mfdataset(hfiles, data_vars="minimal")

            # set timestamp
            base_date = x2d.time.base_date
            year, month, day, hour, tz = base_date.split()
            hour = "{0:05.2f}".format(float(hour)).replace(".", ":")
            tz = "{0:+06.2f}".format(float(tz)).replace(".", "")
            date = "{}-{}-{} {} {}".format(year, month, day, hour, tz)
            # set the start timestamp
            sdate = pd.to_datetime(date, utc=True, format="%Y-%m-%d %H:%M %z")

            # fix fortran/python index
            x2d["SCHISM_hgrid_face_nodes"][:, :3] = x2d["SCHISM_hgrid_face_nodes"].values[:, :3] - 1
            # set time to Datetime
            times = pd.to_datetime(x2d.time.values, unit="s", origin=sdate.tz_convert(None))

            x2d = x2d.assign_coords({"time": ("time", times, x2d.time.attrs)})

            logger.info("get combined 3D netcdf files \n")

            xfiles = glob.glob(os.path.join(path, "outputs/[!out2d_, !hotstart_,]*.nc"))
            xfiles = [x for x in xfiles if not x.endswith("schout_1.nc")]

            if len(xfiles) > 0:
                xfiles.sort()
                # read
                x3d = xr.open_mfdataset(xfiles, data_vars="minimal")
                # set time to Datetime
                x3d = x3d.assign_coords({"time": ("time", times, x3d.time.attrs)})
                x3d.to_netcdf(os.path.join(path, "outputs/schout_2.nc"))

            # save 2D variables to file
            x2d.to_netcdf(os.path.join(path, "outputs/schout_1.nc"))

        else:
            if len(self.misc) == 0:
                logger.info("retrieving index references ... \n")
                self.global2local(**kwargs)
                logger.info("... done \n")

            # Create mesh xarray Dataset
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
                    "SCHISM_hgrid_face_x": (
                        ["nSCHISM_hgrid_face"],
                        gt3.loc[:, "xc"].values,
                    ),
                    "SCHISM_hgrid_face_y": (
                        ["nSCHISM_hgrid_face"],
                        gt3.loc[:, "yc"].values,
                    ),
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
            one = xr.Dataset(
                {
                    "dry_value_flag": (("one"), [iwet_dry]),
                    "SCHISM_hgrid": (("one"), [ihgrid_id]),
                }
            )

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
            hfiles = glob.glob(os.path.join(path, "outputs/schout_*_*.nc"))

            irange_ = [int(x.split("_")[-1].split(".")[0]) for x in hfiles]
            irange_ = np.unique(irange_)

            irange = get_value(self, kwargs, "rlist", irange_)

            self.read_vgrid()  # read grid attributes

            logger.info("write combined netcdf files \n")
            #        total_xdat = []
            for val in tqdm(irange):
                hfiles = glob.glob(os.path.join(path, f"outputs/schout_*_{val}.nc"))
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

                xc.Cs.attrs = {
                    "long_name": "Function C(s) at whole levels",
                    "positive": "up",
                }

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
                xc.time.attrs = {
                    "long_name": "Time",
                    "base_date": base_date,
                    "standard_name": "time",
                }

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

                xc.to_netcdf(os.path.join(path, f"outputs/schout_{val}.nc"))

        logger.info("done with output netCDF files \n")

    def set_obs(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./schism/")
        nspool_sta = get_value(self, kwargs, "nspool_sta", 1)
        tg_database = get_value(self, kwargs, "obs", None)

        if tg_database == None:
            logger.info("get stations using searvey\n")
            geometry = get_value(self, kwargs, "geometry", None)
            if geometry == "global":
                tg = ioc.get_ioc_stations()
            else:
                geo_box = shapely.geometry.box(self.lon_min, self.lat_min, self.lon_max, self.lat_max)
                tg = ioc.get_ioc_stations(region=geo_box)
        else:
            logger.info("get stations from {}\n".format(tg_database))
            tg = pd.read_csv(tg_database)

        tg = tg.reset_index(drop=True)
        ### save in compatible to searvey format
        tg["country"] = tg.country.values.astype("str")  # fix an issue with searvey see #43 therein
        logger.info("save station DataFrame \n")
        sfilename = os.path.join(path, "stations.json")
        tg.to_file(sfilename)
        self.obs = sfilename

        ##### normalize to be used inside pyposeidon
        tgn = normalize_column_names(tg.copy())

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

        logger.info("set in-situ measurements locations \n")

        if not self.mesh.Dataset:
            logger.warning("no mesh available skipping \n")
            return

        gpoints = np.array(
            list(
                zip(
                    self.mesh.Dataset.SCHISM_hgrid_node_x.values,
                    self.mesh.Dataset.SCHISM_hgrid_node_y.values,
                )
            )
        )

        stations = []
        mesh_index = []
        for l in range(tgn.shape[0]):
            plat, plon = tgn.loc[l, ["latitude", "longitude"]]
            cp = closest_node([plon, plat], gpoints)
            mesh_index.append(
                list(
                    zip(
                        self.mesh.Dataset.SCHISM_hgrid_node_x.values,
                        self.mesh.Dataset.SCHISM_hgrid_node_y.values,
                    )
                ).index(tuple(cp))
            )
            stations.append(cp)

        if stations == []:
            logger.warning("no observations available\n")

        # to df
        stations = pd.DataFrame(stations, columns=["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y"])
        stations["z"] = 0
        stations.index += 1
        stations["gindex"] = mesh_index
        try:
            stations["location"] = self.obs.location.values
            stations["provider_id"] = self.obs.ioc_code.values
            stations["provider"] = "ioc"
            stations["longitude"] = self.obs.longitude.values
            stations["latitude"] = self.obs.latitude.values
        except:
            pass

        if coastal_monitoring:
            ## FOR COASTAL MONITORING
            logger.info("set land boundaries as observation points \n")

            # get land boundaries
            coasts = [key for key in self.mesh.Dataset.variables if "land_boundary" in key]
            # get index
            bnodes = []
            for i in range(len(coasts)):
                bnodes = bnodes + list(self.mesh.Dataset[coasts[i]].dropna("index").astype(int).values)
            bnodes = np.unique(bnodes)  # keep unique values
            # get x,y
            xx = self.mesh.Dataset.SCHISM_hgrid_node_x.to_dataframe().loc[bnodes]
            yy = self.mesh.Dataset.SCHISM_hgrid_node_y.to_dataframe().loc[bnodes]
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
        self.params.write(os.path.join(path, "param.nml"), force=True)

        self.stations = stations

        logger.info("write out stations.in file \n")

        # output to file
        with open(os.path.join(path, "station.in"), "w") as f:
            station_flag.to_csv(f, header=None, index=False, sep=" ")
            f.write("{}\n".format(stations.shape[0]))
            stations.loc[:, ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "z"]].to_csv(f, header=None, sep=" ")

    def get_station_sim_data(self, **kwargs):
        path = get_value(self, kwargs, "rpath", "./schism/")

        # locate the station files
        sfiles = glob.glob(os.path.join(path, "outputs/staout_*"))
        sfiles.sort()

        try:
            # get the station flags
            flags = pd.read_csv(os.path.join(path, "station.in"), header=None, nrows=1, delim_whitespace=True).T
            flags.columns = ["flag"]
            flags["variable"] = [
                "elev",
                "air_pressure",
                "windx",
                "windy",
                "T",
                "S",
                "u",
                "v",
                "w",
            ]

            vals = flags[flags.values == 1]  # get the active ones
        except OSError as e:
            if e.errno == errno.EEXIST:
                logger.error("no station.in file present")
            return

        dstamp = kwargs.get("dstamp", self.rdate)

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
            r.index.names = ["time", "node"]

            r.index = r.index.set_levels(r.index.levels[1] - 1, level=1)

            dfs.append(r.to_xarray())

        self.station_sim_data = xr.combine_by_coords(dfs)

    def get_output_data(self, **kwargs):
        dic = self.__dict__

        dic.update(kwargs)

        self.data = data.get_output(**dic)

    def get_station_obs_data(self, **kwargs):
        self.station_obs_data = get_obs_data(self.obs, self.start_date, self.end_date)

    def open_thalassa(self, **kwargs):
        # open a Thalassa instance to visualize the output
        return
