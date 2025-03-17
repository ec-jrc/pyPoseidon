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
from pyposeidon.mesh import r2d
import pyposeidon.model as pm
from pyposeidon.tools import flat_list
import xarray as xr
import glob
import logging
import json

from .. import tools

logger = logging.getLogger(__name__)


def get_output(solver_name: str, **kwargs):
    if solver_name == "schism":
        solver_class = SchismResults
    elif solver_name == "d3d":
        solver_class = D3DResults
    elif solver_name == "telemac":
        solver_class = TelemacResults
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
            data = json.load(f)
            data = pd.json_normalize(data, max_level=0)
            self.info = data.to_dict(orient="records")[0]

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
                    data = json.load(f)
                    data = pd.json_normalize(data, max_level=0)
                    info = data.to_dict(orient="records")[0]

                p = pm.set(**info)

                p.misc = misc

                p.results()

                self.misc = p.misc

                xdat = glob.glob(folder + "/outputs/schout_[!0]*.nc")
                xdat.sort(key=lambda f: int("".join(filter(str.isdigit, f))))

                datai.append(xdat)  # append to list

        if not any(datai):
            logger.warning("no output netcdf files.")
            self.Dataset = None
        else:
            merge = kwargs.get("merge", True)

            if merge:
                datai = flat_list(datai)
                self.Dataset = xr.open_mfdataset(datai, combine="by_coords", data_vars="minimal")

            else:
                self.Dataset = [xr.open_mfdataset(x, combine="by_coords", data_vars="minimal") for x in datai]


class TelemacResults:
    def __init__(self, **kwargs):
        """
        this class has been copied on the Schism class above
        although there are a few subtilities that TELEMAC introduced :
         1. there are 1D AND 2D results, they need to be addressed separately
         2. selafin are readable directly in xarray, so we can skip the conversion step
        """
        from pyposeidon.telemac import extract_t_elev_2D

        rpath = kwargs.get("rpath", "./telemac/")
        res_type = kwargs.get("result_type", "2D")
        convert = kwargs.get("convert_results", True)
        extract_TS = kwargs.get("extract_TS", False)
        id_str = kwargs.get("id_str", "ioc_code")
        max_dist = kwargs.get("max_dist", 1000)

        if res_type not in ["1D", "2D"]:
            raise ValueError("results_type needs to be '1D' or '2D'!")
        if res_type == "1D":
            out_default = "stations.zarr"
        else:
            out_default = "out_2D.zarr"

        folders = kwargs.get("folders", None)

        if folders:
            self.folders = folders
        else:
            self.folders = [rpath]

        datai = []

        module = kwargs.get("module", "telemac2d")

        misc = kwargs.get("misc", {})

        for folder in self.folders:
            logger.info(" Combining output for folder {}\n".format(folder))

            with open(folder + "/" + module + "_model.json", "r") as f:
                data = json.load(f)
                data = pd.json_normalize(data, max_level=0)
                info = data.to_dict(orient="records")[0]

            p = pm.set(**info)
            p.misc = misc

            # read from output file
            if convert:
                xdat = glob.glob(folder + "/outputs/" + out_default)
                model_xstr, model_ystr = "longitude", "latitude"
                if module == "telemac2d":
                    var = "elev"
                elif module == "telemac3d":
                    var = "elev"
                elif module == "tomawac":
                    var = "hm0"
                else:
                    raise ValueError(f"module {module} not supported!")

                if len(xdat) > 0:
                    datai.append(xdat)  # append to list

                else:  # run merge output
                    p.results(
                        filename="stations.zarr",
                        filename2d="out_2D.zarr",
                        remove_zarr=False,  # remove zarr files after tarballing
                    )

                    self.misc = p.misc
                    xdat = glob.glob(folder + "/outputs/" + out_default)
                    datai.append(xdat)  # append to list
            else:  # read from selafin file
                xdat = glob.glob(folder + f"/results_{res_type}.slf")
                datai.append(xdat)  # append to list
                model_xstr, model_ystr = "x", "y"
                if module == "telemac2d":
                    var = "S"
                elif module == "telemac3d":
                    var = "Z"
                elif module == "tomawac":
                    var = "WH"
                else:
                    raise ValueError(f"module {module} not supported!")

        if not any(datai):
            logger.warning("no output files found")
            self.Dataset = None
        else:
            merge = kwargs.get("merge", True)

            if merge:
                datai = flat_list(datai)
                self.Dataset = xr.open_mfdataset(datai, data_vars="minimal")
            else:
                self.Dataset = [xr.open_mfdataset(x, data_vars="minimal") for x in datai]

        # export parquet time series
        if extract_TS:
            if "stations" in p.__dict__:
               file = p.stations
            elif "stations.csv" in os.listdir(rpath):
                file = os.path.join(rpath, "stations.csv")
            elif "obs" in self.__dict__:
                p.set_obs()
                file = p.stations
            else:
                raise ValueError("no stations file info found")
            stations = pd.read_csv(file)
            logger.info("extracting parquet files from TELEMAC Selafin output \n")
            for i_s, id_ in enumerate(stations[id_str]):
                s = stations[stations[id_str] == id_]
                mod, mlon, mlat = extract_t_elev_2D(
                    self.Dataset,
                    s.lon.values[0],
                    s.lat.values[0],
                    var,
                    model_xstr,
                    model_ystr,
                    max_dist=max_dist,
                )
                mod.to_frame().to_parquet(f"{rpath}{id_}.parquet")
