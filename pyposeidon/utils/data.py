"""
Data analysis module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import pandas as pd
import os
from pyposeidon.utils.vals import obs
from pyposeidon.mesh import r2d
import pyposeidon.model as pm
from pyposeidon.tools import flat_list
from pyposeidon.utils.get_value import get_value
import datetime
import xarray as xr
import glob
import sys
import logging

from .. import tools

logger = logging.getLogger(__name__)


def get_output(solver_name: str, **kwargs):
    if solver_name == "schism":
        solver_class = SchismResults
    elif solver_name == "d3d":
        solver_class = D3DResults
    else:
        raise ValueError(f"Unknown solver_name: {solver_name}")
    instance = solver_class(**kwargs)
    return instance


class D3DResults:
    def __init__(self, **kwargs):

        rpath = kwargs.get("rpath", "./d3d/")

        folders = kwargs.get(
            "folders", None
        )  # [os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]

        if folders:
            self.folders = folders
        else:
            self.folders = [rpath]

        # check if many tags present
        ifiles = glob.glob(self.folders[0] + "/*_model.json")

        if len(ifiles) > 1:
            # ---------------------------------------------------------------------
            logger.warning("more than one configuration, specify tag argument \n")
        # ---------------------------------------------------------------------

        tag = kwargs.get("tag", None)

        if tag:
            ifile = self.folders[0] + "/" + tag + "_model.json"
        else:
            ifile = ifiles[0]

        # ---------------------------------------------------------------------
        logger.info("reading data based on {} \n".format(ifile))
        # ---------------------------------------------------------------------

        with open(ifile, "rb") as f:
            info = pd.read_json(f, lines=True).T
            info[info.isnull().values] = None
            self.info = info.to_dict()[0]

        grid = r2d.read_file(self.folders[0] + "/" + self.info["tag"] + ".grd")

        deb = np.loadtxt(self.folders[0] + "/" + self.info["tag"] + ".dep")

        # create mask

        d = deb[1:-1, 1:-1]
        self.w = d == -999.0

        b = deb[:-1, :-1]
        b[b == -999.0] = np.nan

        self.dem = xr.Dataset(
            {"bathymetry": (["latitude", "longitude"], -b)},
            coords={
                "longitude": ("longitude", grid.lons[0, :].values),
                "latitude": ("latitude", grid.lats[:, 0].values),
            },
        )

        self.grid = grid

        # READ DATA

        nfiles = [folder + "/" + "trim-" + self.info["tag"] + ".nc" for folder in self.folders]

        ds = xr.open_mfdataset(nfiles, combine="by_coords", data_vars="minimal")

        self.Dataset = ds

        # clean duplicates
        self.Dataset = self.Dataset.sel(time=~self.Dataset.indexes["time"].duplicated())

        dic = self.info.copy()  # start with x's keys and values
        dic.update(kwargs)  # modifies z with y's keys and values & returns None

        if "sa_date" not in dic.keys():
            dic.update({"sa_date": self.Dataset.time.values[0]})

        if "se_date" not in dic.keys():
            dic.update({"se_date": self.Dataset.time.values[-1]})

        self.obs = obs(**dic)

    def frames(self, var, **kwargs):

        X, Y = self.Dataset.XZ.values[1:-1, 1:-1], self.Dataset.YZ.values[1:-1, 1:-1]
        xh = np.ma.masked_array(X.T, self.w)  # mask land
        yh = np.ma.masked_array(Y.T, self.w)

        if len(var) == 1:

            var = self.Dataset[var[0]].transpose(
                self.Dataset[var[0]].dims[0],
                self.Dataset[var[0]].dims[2],
                self.Dataset[var[0]].dims[1],
                transpose_coords=True,
            )[:, 1:-1, 1:-1]
            ww = np.broadcast_to(self.w == True, var.shape)
            v = np.ma.masked_array(var, ww)
            return contour(xh, yh, v, self.Dataset.time.values, **kwargs)

        elif len(var) == 2:

            a0 = self.Dataset[var[0]].squeeze()
            var0 = a0.transpose(a0.dims[0], a0.dims[2], a0.dims[1])[:, 1:-1, 1:-1]
            a1 = self.Dataset[var[1]].squeeze()
            var1 = a1.transpose(a1.dims[0], a1.dims[2], a1.dims[1])[:, 1:-1, 1:-1]
            ww = np.broadcast_to(self.w == True, var0.shape)

            v0 = np.ma.masked_array(var0, ww)
            v1 = np.ma.masked_array(var1, ww)
            return quiver(xh, yh, v0, v1, self.Dataset.time.values, **kwargs)


class SchismResults:
    def __init__(self, **kwargs):

        rpath = kwargs.get("rpath", "./schism/")

        folders = kwargs.get(
            "folders", None
        )  # [os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]

        if folders:
            self.folders = folders
        else:
            self.folders = [rpath]

        datai = []

        tag = kwargs.get("tag", "schism")

        misc = kwargs.get("misc", {})

        for folder in self.folders:

            logger.info(" Combining output for folder {}\n".format(folder))

            xdat = glob.glob(folder + "/outputs/schout_[!0]*.nc")
            xdat.sort(key=lambda f: int("".join(filter(str.isdigit, f))))

            if len(xdat) > 0:
                datai.append(xdat)  # append to list

            else:  # run merge output

                with open(folder + "/" + tag + "_model.json", "r") as f:
                    info = pd.read_json(f, lines=True).T
                    info[info.isnull().values] = None
                    info = info.to_dict()[0]

                p = pm.set(**info)

                p.misc = misc

                p.results()

                self.misc = p.misc

                xdat = glob.glob(folder + "/outputs/schout_[!0]*.nc")
                xdat.sort(key=lambda f: int("".join(filter(str.isdigit, f))))

                datai.append(xdat)  # append to list

        merge = kwargs.get("merge", True)

        if merge:

            datai = flat_list(datai)
            self.Dataset = xr.open_mfdataset(datai, combine="by_coords", data_vars="minimal")

            with open(self.folders[-1] + "/" + tag + "_model.json", "r") as f:
                info = pd.read_json(f, lines=True).T
                info[info.isnull().values] = None
                info = info.to_dict()[0]

            p = pm.set(**info)

            if hasattr(p, "stations"):

                logger.info("Retrieve station timeseries\n")

                dstamp = kwargs.get("dstamp", info["date"])

                p.get_station_data(dstamp=dstamp)
                self.time_series = p.time_series

        else:
            self.Dataset = [xr.open_mfdataset(x, combine="by_coords", data_vars="minimal") for x in datai]

            ts = []

            for folder in self.folders:

                p = pm.read_model(folder + "/{}_model.json".format(tag))  # read model

                if hasattr(p, "stations"):

                    logger.info("Retrieve station timeseries\n")

                    dstamp = kwargs.get("dstamp", p.date)

                    p.get_station_data(dstamp=dstamp)
                    ts.append(p.time_series)

            self.time_series = ts
