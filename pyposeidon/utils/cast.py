"""
Simulation management module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.
import pyposeidon
import pyposeidon.model as pm
import pyposeidon.mesh as pmesh
from pyposeidon.utils.get_value import get_value


import numpy as np
import errno
import datetime
import sys
import os, errno
from shutil import copy2
import glob
import pandas as pd
import pathlib
import json

import logging

logger = logging.getLogger(__name__)


def set(solver_name: str, **kwargs):
    if solver_name == "d3d":
        instance = D3DCast(**kwargs)
    elif solver_name == "schism":
        instance = SchismCast(**kwargs)
    else:
        raise ValueError(f"Don't know how to handle solver: {solver_name}")
    return instance


def copy_files(rpath: str, ppath: str, filenames: list[str]) -> None:
    for filename in filenames:
        src = os.path.join(ppath, filename)
        dst = os.path.join(rpath, filename)
        if os.path.exists(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            copy2(src, dst)
            logger.debug("copied src -> dst: %s -> %s", src, dst)


def symlink_files(rpath: str, ppath: str, filenames: list[str]) -> None:
    for filename in filenames:
        src = os.path.join(ppath, filename)
        dst = os.path.join(rpath, filename)
        if os.path.exists(src):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)
            logger.debug("symlinked src -> dst: %s -> %s", src, dst)


class D3DCast:
    def __init__(self, **kwargs):
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def run(self, **kwargs):
        if isinstance(self.model, str):
            self.model = pyposeidon.model.read(self.model)

        for attr, value in self.model.__dict__.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

        execute = get_value(self, kwargs, "execute", False)

        pwd = os.getcwd()

        files = [
            self.tag + "_hydro.xml",
            self.tag + ".enc",
            self.tag + ".obs",
            self.tag + ".bnd",
            self.tag + ".bca",
            "run_flow2d3d.sh",
        ]
        files_sym = [self.tag + ".grd", self.tag + ".dep"]

        self.origin = self.model.rpath
        self.rdate = self.model.rdate

        if not os.path.exists(self.origin):
            sys.stdout.write("Initial folder not present {}\n".format(self.origin))
            sys.exit(1)

        ppath = self.ppath

        cf = [glob.glob(ppath + "/" + e) for e in files]
        cfiles = [item.split("/")[-1] for sublist in cf for item in sublist]

        # create the folder/run path

        rpath = self.cpath

        if not os.path.exists(rpath):
            os.makedirs(rpath)

        copy2(ppath + self.tag + "_model.json", rpath)  # copy the info file

        # load model
        with open(rpath + self.tag + "_model.json", "rb") as f:
            data = json.load(f)
            data = pd.json_normalize(data, max_level=0)
            info = data.to_dict(orient="records")[0]

        try:
            args = set(kwargs.keys()).intersection(info.keys())  # modify dic with kwargs
            for attr in list(args):
                info[attr] = kwargs[attr]
        except:
            pass

        # update the properties
        info["rdate"] = self.rdate
        info["start_date"] = self.start_date
        info["time_frame"] = self.time_frame
        info["meteo_source"] = self.meteo
        info["rpath"] = rpath
        if self.restart_step:
            info["restart_step"] = self.restart_step

        m = pm.set(**info)

        # copy/link necessary files
        logger.debug("copy necessary files")

        for filename in cfiles:
            ipath = glob.glob(ppath + filename)
            if ipath:
                try:
                    copy2(ppath + filename, rpath + filename)
                except:
                    dir_name, file_name = os.path.split(filename)
                    if not os.path.exists(rpath + dir_name):
                        os.makedirs(rpath + dir_name)
                    copy2(ppath + filename, rpath + filename)
        logger.debug(".. done")

        # symlink the big files
        logger.debug("symlink model files")
        for filename in files_sym:
            ipath = glob.glob(os.path.join(self.origin, filename))
            if ipath:
                try:
                    os.symlink(pathlib.Path(ipath[0]).resolve(strict=True), rpath + filename)
                except OSError as e:
                    if e.errno == errno.EEXIST:
                        logger.warning("Restart link present\n")
                        logger.warning("overwriting\n")
                        os.remove(rpath + filename)
                        os.symlink(
                            pathlib.Path(ipath[0]).resolve(strict=True),
                            rpath + filename,
                        )
        logger.debug(".. done")

        copy2(ppath + m.tag + ".mdf", rpath)  # copy the mdf file

        # copy restart file

        inresfile = "tri-rst." + m.tag + "." + datetime.datetime.strftime(self.rdate, "%Y%m%d.%H%M%M")

        outresfile = "restart." + datetime.datetime.strftime(self.rdate, "%Y%m%d.%H%M%M")

        #  copy2(ppath+inresfile,rpath+'tri-rst.'+outresfile)
        try:
            os.symlink(
                pathlib.Path(ppath + "/" + inresfile).resolve(strict=True),
                rpath + "tri-rst." + outresfile,
            )
            logger.debug("symlink {} to {}".format(ppath + "/" + inresfile, rpath + "tri-rst." + outresfile))
        except OSError as e:
            if e.errno == errno.EEXIST:
                logger.warning("Restart symlink present\n")
                logger.warning("overwriting\n")
                os.remove(rpath + "tri-rst." + outresfile)
                os.symlink(
                    pathlib.Path(ppath + "/" + inresfile).resolve(strict=True),
                    rpath + "tri-rst." + outresfile,
                )
            else:
                raise e

        # get new meteo

        logger.info("process meteo\n")

        flag = get_value(self, kwargs, "update", ["meteo"])

        check = [os.path.exists(rpath + f) for f in ["u.amu", "v.amv", "p.amp"]]

        if (np.any(check) == False) or ("meteo" in flag):
            m.force()
            m.to_force(m.meteo.Dataset, vars=["msl", "u10", "v10"], rpath=rpath)  # write u,v,p files

        else:
            logger.info("meteo files present\n")

        # modify mdf file
        m.config(
            config_file=ppath + m.tag + ".mdf",
            config={"Restid": outresfile},
            output=True,
        )

        m.config_file = rpath + m.tag + ".mdf"

        os.chdir(rpath)
        m.save()

        if execute:
            m.run()

        # cleanup
        os.remove(rpath + "tri-rst." + outresfile)

        logger.info("done for date :" + datetime.datetime.strftime(self.rdate, "%Y%m%d.%H"))

        os.chdir(pwd)


class SchismCast:
    files = [
        "launchSchism.sh",
        "sflux/sflux_inputs.txt",
        "outputs/flux.out",
    ]

    model_files = [
        "bctides.in",
        "hgrid.gr3",
        "hgrid.ll",
        "manning.gr3",
        "vgrid.in",
        "drag.gr3",
        "rough.gr3",
        "station.in",
        "stations.json",
        "windrot_geo2proj.gr3",
    ]

    station_files = [
        "outputs/staout_1",
        "outputs/staout_2",
        "outputs/staout_3",
        "outputs/staout_4",
        "outputs/staout_5",
        "outputs/staout_6",
        "outputs/staout_7",
        "outputs/staout_8",
        "outputs/staout_9",
    ]

    def __init__(self, **kwargs):
        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def run(self, **kwargs):
        if isinstance(self.model, str):
            self.model = pyposeidon.model.read(self.model)

        for attr, value in self.model.__dict__.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

        execute = get_value(self, kwargs, "execute", True)

        copy = get_value(self, kwargs, "copy", False)

        pwd = os.getcwd()

        self.origin = self.model.rpath
        self.rdate = self.model.rdate
        ppath = self.ppath

        # control
        if not isinstance(self.rdate, pd.Timestamp):
            self.rdate = pd.to_datetime(self.rdate)

        if not os.path.exists(self.origin):
            sys.stdout.write(f"Initial folder not present {self.origin}\n")
            sys.exit(1)

        # create the new folder/run path
        rpath = self.cpath

        if not os.path.exists(rpath):
            os.makedirs(rpath)

        model_definition_filename = f"{self.tag}_model.json"
        copy2(os.path.join(ppath, model_definition_filename), rpath)  # copy the info file

        # load model
        with open(os.path.join(rpath, model_definition_filename), "rb") as f:
            data = json.load(f)
            data = pd.json_normalize(data, max_level=0)
            info = data.to_dict(orient="records")[0]

        try:
            args = set(kwargs.keys()).intersection(info.keys())  # modify dic with kwargs
            for attr in list(args):
                info[attr] = kwargs[attr]
        except:
            pass

        # add optional additional kwargs
        for attr in kwargs.keys():
            if attr not in info.keys():
                info[attr] = kwargs[attr]

        info["config_file"] = os.path.join(ppath, "param.nml")

        # update the properties

        info["rdate"] = self.rdate
        info["start_date"] = self.sdate
        info["time_frame"] = self.time_frame
        info["end_date"] = self.sdate + pd.to_timedelta(self.time_frame)
        info["meteo_source"] = self.meteo
        info["rpath"] = rpath

        m = pm.set(**info)

        # Mesh
        gfile = glob.glob(os.path.join(ppath, "hgrid.gr3"))
        if gfile:
            info["mesh_file"] = gfile[0]
            self.mesh_file = gfile[0]
            info["mesh_generator"] = None
            self.mesh_generator = None

        m.mesh = pmesh.set(type="tri2d", **info)

        # get lat/lon from file
        if hasattr(self, "mesh_file"):
            info.update({"lon_min": m.mesh.Dataset.SCHISM_hgrid_node_x.values.min()})
            info.update({"lon_max": m.mesh.Dataset.SCHISM_hgrid_node_x.values.max()})
            info.update({"lat_min": m.mesh.Dataset.SCHISM_hgrid_node_y.values.min()})
            info.update({"lat_max": m.mesh.Dataset.SCHISM_hgrid_node_y.values.max()})

        # copy/link necessary files
        logger.debug("Copy necessary + station files")
        copy_files(rpath=rpath, ppath=ppath, filenames=self.files + self.station_files)
        logger.debug("Copy model files")
        if copy:
            copy_files(rpath=rpath, ppath=ppath, filenames=self.model_files)
        else:
            symlink_files(rpath=rpath, ppath=ppath, filenames=self.model_files)

        logger.debug(".. done")

        # create restart file
        logger.debug("create restart file")

        # check for combine hotstart
        hotout = int((self.sdate - self.rdate).total_seconds() / info["params"]["core"]["dt"])
        logger.debug("hotout_it = {}".format(hotout))

        # link restart file
        inresfile = os.path.join(ppath, f"outputs/hotstart_it={hotout}.nc")
        outresfile = os.path.join(rpath, "hotstart.nc")

        logger.debug("hotstart_file: %s", inresfile)
        if not os.path.exists(inresfile):
            logger.info("Generating hotstart file.\n")
            # load model model from ppath
            with open(os.path.join(ppath, self.tag + "_model.json"), "rb") as f:
                data = json.load(f)
                data = pd.json_normalize(data, max_level=0)
                ph = data.to_dict(orient="records")[0]
            p = pm.set(**ph)
            p.hotstart(it=hotout)
        else:
            logger.info("Hotstart file already existing. Skipping creation.\n")

        if copy:
            logger.info("Copying: %s -> %s", inresfile, outresfile)
            copy2(inresfile, outresfile)
        else:
            logger.info("Symlinking`: %s -> %s", inresfile, outresfile)
            try:
                os.symlink(
                    pathlib.Path(os.path.join(ppath, inresfile)).resolve(strict=True), os.path.join(rpath, outresfile)
                )
            except OSError as e:
                if e.errno == errno.EEXIST:
                    logger.warning("Restart link present\n")
                    logger.warning("overwriting\n")
                    os.remove(os.path.join(rpath, outresfile))
                    os.symlink(
                        pathlib.Path(os.path.join(ppath, inresfile)).resolve(strict=True),
                        os.path.join(rpath, outresfile),
                    )
                else:
                    raise e
        # get new meteo

        logger.info("process meteo\n")

        flag = get_value(self, kwargs, "update", [])

        check = [os.path.exists(os.path.join(rpath, "sflux", f)) for f in ["sflux_air_1.0001.nc"]]

        if (np.any(check) == False) or ("meteo" in flag):
            m.force(**info)
            if hasattr(self, "meteo_split_by"):
                times, datasets = zip(*m.meteo.Dataset.groupby("time.{}".format(self.meteo_split_by)))
                mpaths = ["sflux_air_1.{:04d}.nc".format(t + 1) for t in np.arange(len(times))]
                for das, mpath in list(zip(datasets, mpaths)):
                    m.to_force(
                        das,
                        vars=["msl", "u10", "v10"],
                        rpath=rpath,
                        filename=mpath,
                        date=self.rdate,
                    )
            else:
                m.to_force(
                    m.meteo.Dataset,
                    vars=["msl", "u10", "v10"],
                    rpath=rpath,
                    date=self.rdate,
                )

        else:
            logger.warning("meteo files present\n")

        # modify param file
        rnday_new = (self.sdate - self.rdate).total_seconds() / (3600 * 24.0) + pd.to_timedelta(
            self.time_frame
        ).total_seconds() / (3600 * 24.0)
        hotout_write = int(rnday_new * 24 * 3600 / info["params"]["core"]["dt"])
        info["parameters"].update(
            {
                "ihot": 2,
                "rnday": rnday_new,
                "start_hour": self.rdate.hour,
                "start_day": self.rdate.day,
                "start_month": self.rdate.month,
                "start_year": self.rdate.year,
            }
        )

        m.config(output=True, **info)  # save param.nml

        m.config_file = os.path.join(rpath, "param.nml")

        m.save()

        if execute:
            m.run()

        logger.info("done for date : %s", self.sdate.strftime("%Y%m%d.%H"))

        os.chdir(pwd)
