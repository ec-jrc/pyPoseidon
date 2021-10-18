"""
Data analysis module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import os
from pyposeidon.utils.vals import obs
from pyposeidon.grid import *
import pyposeidon.model as pm
import datetime
from pyposeidon.utils.get_value import get_value
import xarray as xr
import glob
import sys
import logging

logger = logging.getLogger("pyposeidon")


def data(**kwargs):

    solver = kwargs.get("solver", None)
    if solver == "d3d":
        return d3d(**kwargs)
    elif solver == "schism":
        return schism(**kwargs)
    else:
        logger.error("solver is not defined, exiting \n")
        sys.exit(1)


class d3d:
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
                "longitude": ("longitude", grid.lons[0, :]),
                "latitude": ("latitude", grid.lats[:, 0]),
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


class schism:
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

            logger.info("Combining output for folder {}\n".format(folder))

            xdat = glob.glob(folder + "/outputs/schout_[!0]*.nc")
            xdat.sort(key=lambda f: int("".join(filter(str.isdigit, f))))

            if len(xdat) > 0:
                datai.extend(xdat)  # append to list

            else:  # run merge output

                with open(folder + "/" + tag + "_model.json", "r") as f:
                    info = pd.read_json(f, lines=True).T
                    info[info.isnull().values] = None
                    info = info.to_dict()[0]

                p = pm.model(**info)

                p.misc = misc

                p.results()

                self.misc = p.misc

                xdat = glob.glob(folder + "/outputs/schout_[!0]*.nc")
                xdat.sort(key=lambda f: int("".join(filter(str.isdigit, f))))

                datai.extend(xdat)  # append to list

        self.Dataset = xr.open_mfdataset(datai, combine="by_coords", data_vars="minimal")

        with open(self.folders[-1] + "/" + tag + "_model.json", "r") as f:
            info = pd.read_json(f, lines=True).T
            info[info.isnull().values] = None
            info = info.to_dict()[0]

        p = pm.model(**info)

        logger.info("Retrieve station timeseries if any\n")

        dstamp = kwargs.get("dstamp", info["date"])
        try:
            p.get_station_data(dstamp=dstamp)
            self.time_series = p.time_series

            dic = {}

            try:

                with open(self.folders[0] + "/" + tag + "_model.json", "r") as f:
                    info = pd.read_json(f, lines=True).T
                    info[info.isnull().values] = None
                    self.info = info.to_dict()[0]

                dic = self.info.copy()  # start with x's keys and values
                dic.update(kwargs)  # modifies z with y's keys and values & returns None

            except:
                pass

            if "sa_date" not in dic.keys():
                dic.update({"sa_date": self.time_series.time.values[0]})

            if "se_date" not in dic.keys():
                dic.update({"se_date": self.Dataset.time.values[-1]})

            if "lon_min" not in dic.keys():
                dic.update({"lon_min": self.Dataset.SCHISM_hgrid_node_x.values.min()})

            if "lon_max" not in dic.keys():
                dic.update({"lon_max": self.Dataset.SCHISM_hgrid_node_x.values.max()})

            if "lat_min" not in dic.keys():
                dic.update({"lat_min": self.Dataset.SCHISM_hgrid_node_y.values.min()})

            if "lat_max" not in dic.keys():
                dic.update({"lat_max": self.Dataset.SCHISM_hgrid_node_y.values.max()})

            logger.info("Retrieve observations info\n")

            self.obs = obs(**dic)

            ret = kwargs.get("online", False)

            if ret is True:
                logger.info("collect observational data")
                tgs = self.obs.locations.loc[self.obs.locations.Group == "TD UNESCO"]

                dic = {}
                for i in tgs.index:
                    #            print(i, tgs.loc[i].Name.strip())
                    while True:
                        p = self.obs.iloc(i)
                        if p is not None:
                            if p.shape[0] > 1:
                                p = p.dropna()
                            break
                    dic.update({tgs.loc[i].Name.strip(): p})

                tg = pd.concat(dic, axis=0, sort=True)

                try:
                    tg = tg.drop("TimeUTC", axis=1)
                except:
                    pass

                tg.to_csv(self.folders[0] + "/" + "obs.csv")

                self.obs.dataframe = tg

        except:
            logger.info("no station data loaded")
