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

# from pyposeidon.utils import data
import subprocess
import logging

logger = logging.getLogger("pyposeidon")


def set(solver=None, **kwargs):
    if solver == "d3d":
        return dcast(**kwargs)
    elif solver == "schism":
        return scast(**kwargs)


class dcast:
    def __init__(self, **kwargs):

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def run(self, **kwargs):

        if isinstance(self.model, str):
            self.model = pyposeidon.read_model(self.model)

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
        self.date0 = self.model.date

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
            info = pd.read_json(f, lines=True).T
            info[info.isnull().values] = None
            info = info.to_dict()[0]

        args = set(kwargs.keys()).intersection(info.keys())  # modify dic with kwargs
        for attr in list(args):
            info[attr] = kwargs[attr]

        # update the properties
        info["date"] = self.date
        info["start_date"] = self.date
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
            ipath = glob.glob(self.origin + filename)
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

        inresfile = "tri-rst." + m.tag + "." + datetime.datetime.strftime(self.date, "%Y%m%d.%H%M%M")

        outresfile = "restart." + datetime.datetime.strftime(self.date, "%Y%m%d.%H%M%M")

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

        logger.info("done for date :" + datetime.datetime.strftime(self.date, "%Y%m%d.%H"))

        os.chdir(pwd)


class scast:
    def __init__(self, **kwargs):

        for attr, value in kwargs.items():
            setattr(self, attr, value)

    def run(self, **kwargs):

        if isinstance(self.model, str):
            self.model = pyposeidon.read_model(self.model)

        for attr, value in self.model.__dict__.items():
            if not hasattr(self, attr):
                setattr(self, attr, value)

        execute = get_value(self, kwargs, "execute", True)

        pwd = os.getcwd()

        files = [
            "bctides.in",
            "launchSchism.sh",
            "/sflux/sflux_inputs.txt",
            "/outputs/flux.out",
        ]
        files_sym = [
            "hgrid.gr3",
            "hgrid.ll",
            "manning.gr3",
            "vgrid.in",
            "drag.gr3",
            "rough.gr3",
            "station.in",
            "windrot_geo2proj.gr3",
        ]
        station_files = [
            "/outputs/staout_1",
            "/outputs/staout_2",
            "/outputs/staout_3",
            "/outputs/staout_4",
            "/outputs/staout_5",
            "/outputs/staout_6",
            "/outputs/staout_7",
            "/outputs/staout_8",
            "/outputs/staout_9",
        ]

        self.origin = self.model.rpath
        self.date0 = self.model.date

        if not os.path.exists(self.origin):
            sys.stdout.write("Initial folder not present {}\n".format(self.origin))
            sys.exit(1)

        ppath = self.ppath
        # create the new folder/run path
        rpath = self.cpath

        if not os.path.exists(rpath):
            os.makedirs(rpath)

        tag = kwargs.get("tag", "schism")
        copy2(ppath + self.tag + "_model.json", rpath)  # copy the info file

        # load model
        with open(rpath + self.tag + "_model.json", "rb") as f:
            info = pd.read_json(f, lines=True).T
            info[info.isnull().values] = None
            info = info.to_dict()[0]

        try:
            args = set(kwargs.keys()).intersection(info.keys())  # modify dic with kwargs
            for attr in list(args):
                info[attr] = kwargs[attr]
        except:
            pass

        info["config_file"] = ppath + "param.nml"

        # update the properties

        info["date"] = self.date0
        info["start_date"] = self.sdate
        info["time_frame"] = self.time_frame
        info["end_date"] = self.sdate + pd.to_timedelta(self.time_frame)
        info["meteo_source"] = self.meteo
        info["rpath"] = rpath

        m = pm.set(**info)

        # Mesh
        gfile = glob.glob(ppath + "hgrid.gr3")
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
        logger.debug("copy necessary files")

        for filename in files:
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

        # copy the station files
        logger.debug("copy station files")
        for filename in station_files:
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
            ipath = glob.glob(self.origin + filename)
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

        # create restart file
        logger.debug("create restart file")

        # check for combine hotstart
        hotout = int((self.sdate - self.date0).total_seconds() / info["params"]["core"]["dt"])
        logger.debug("hotout_it = {}".format(hotout))

        resfile = glob.glob(ppath + "/outputs/hotstart_it={}.nc".format(hotout))
        if not resfile:
            # load model model from ppath
            with open(ppath + self.tag + "_model.json", "rb") as f:
                ph = pd.read_json(f, lines=True).T
                ph[ph.isnull().values] = None
                ph = ph.to_dict()[0]
            p = pm.set(**ph)
            p.hotstart(it=hotout)

        # link restart file
        inresfile = "/outputs/hotstart_it={}.nc".format(hotout)
        outresfile = "/hotstart.nc"

        logger.info("set restart\n")

        try:
            os.symlink(pathlib.Path(ppath + inresfile).resolve(strict=True), rpath + outresfile)
        except OSError as e:
            if e.errno == errno.EEXIST:
                logger.warning("Restart link present\n")
                logger.warning("overwriting\n")
                os.remove(rpath + outresfile)
                os.symlink(
                    pathlib.Path(ppath + inresfile).resolve(strict=True),
                    rpath + outresfile,
                )
            else:
                raise e

        # get new meteo

        logger.info("process meteo\n")

        flag = get_value(self, kwargs, "update", [])

        check = [os.path.exists(rpath + "sflux/" + f) for f in ["sflux_air_1.0001.nc"]]

        if (np.any(check) == False) or ("meteo" in flag):

            m.force(**info)
            if hasattr(self, "meteo_split_by"):
                times, datasets = zip(*m.meteo.Dataset.groupby("time.{}".format(self.meteo_split_by)))
                mpaths = ["sflux/sflux_air_1.{:04d}.nc".format(t + 1) for t in np.arange(len(times))]
                for das, mpath in list(zip(datasets, mpaths)):
                    m.to_force(
                        das,
                        vars=["msl", "u10", "v10"],
                        rpath=rpath,
                        filename=mpath,
                        date=self.date0,
                    )
            else:
                m.to_force(
                    m.meteo.Dataset,
                    vars=["msl", "u10", "v10"],
                    rpath=rpath,
                    date=self.date0,
                )

        else:
            logger.warning("meteo files present\n")

        # modify param file
        rnday_new = (self.sdate - self.date0).total_seconds() / (3600 * 24.0) + pd.to_timedelta(
            self.time_frame
        ).total_seconds() / (3600 * 24.0)
        hotout_write = int(rnday_new * 24 * 3600 / info["params"]["core"]["dt"])
        info["parameters"].update(
            {
                "ihot": 2,
                "rnday": rnday_new,
                "start_hour": self.date0.hour,
                "start_day": self.date0.day,
                "start_month": self.date0.month,
                "start_year": self.date0.year,
            }
        )

        m.config(output=True, **info)  # save param.nml

        m.config_file = rpath + "param.nml"

        m.save()

        if execute:
            m.run()

        logger.info("done for date :" + self.sdate.strftime("%Y%m%d.%H"))

        os.chdir(pwd)
