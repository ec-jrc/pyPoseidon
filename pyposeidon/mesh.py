"""
Mesh module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import datetime
import xarray as xr
import pandas as pd
import sys
from pyposeidon import mjigsaw
from pyposeidon import mgmsh
import logging
import f90nml
import os
import subprocess
from pyposeidon.utils.verify import *
import pyposeidon.boundary as pb

from . import tools

logger = logging.getLogger("pyposeidon")

DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/misc/"


def set(type=None, **kwargs):
    if type == "r2d":
        return r2d(**kwargs)
    elif type == "tri2d":
        return tri2d(**kwargs)


class r2d:
    """Regular 2d grid for d3d"""

    def __init__(self, **kwargs):

        mesh_file = kwargs.get("mesh_file", None)

        if mesh_file:

            self.Dataset = self.read_file(mesh_file)

        else:

            geometry = kwargs.get("geometry", None)

            try:
                lon_min = geometry["lon_min"]
                lon_max = geometry["lon_max"]
                lat_min = geometry["lat_min"]
                lat_max = geometry["lat_max"]
            except:
                logger.error("geometry not set properly\n")
                sys.exit(1)

            resolution = kwargs.get("resolution", None)

            ni = int(round((lon_max - lon_min) / resolution))  # these are cell numbers
            nj = int(round((lat_max - lat_min) / resolution))

            lon_max = lon_min + ni * resolution  # adjust max lon to much the grid
            lat_max = lat_min + nj * resolution

            # set the grid
            x = np.linspace(lon_min, lon_max, ni)
            y = np.linspace(lat_min, lat_max, nj)
            gx, gy = np.meshgrid(x, y)

            attrs = kwargs.get(
                "attrs",
                {
                    "Coordinate System": "Spherical",
                    "alfori": 0.0,
                    "xori": 0.0,
                    "yori": 0.0,
                },
            )

            g = xr.Dataset(
                {"lons": (["y", "x"], gx), "lats": (["y", "x"], gy)},
                coords={"x": ("x", gx[0, :]), "y": ("y", gy[:, 0])},
            )

            g.attrs = attrs

            self.Dataset = g

    @staticmethod
    def read_file(filename, **kwargs):

        logger.info("read grid file {}".format(filename))

        header = pd.read_csv(filename, nrows=3, header=None, comment="*")
        cs = header.loc[0, 0].split("=")[1].strip()
        ni, nj = header.loc[1, 0].split(" ")
        ni, nj = int(ni), int(nj)
        alfori, xori, yori = header.loc[2, 0].split(" ")

        d = pd.read_csv(
            filename,
            header=2,
            comment="*",
            delim_whitespace=True,
            engine="python",
            na_values="ETA=",
        )
        d = d.reset_index()
        data = d.values[~np.isnan(d.values)]
        data = np.array(data)
        data = data.reshape(2, nj, ni + 1)  # including the row index
        # clean up the row index
        data = data[:, :, 1:]

        lons = data[0, :, :]
        lats = data[1, :, :]

        g = xr.Dataset(
            {"lons": (["y", "x"], lons), "lats": (["y", "x"], lats)},
            coords={"x": ("x", lons[0, :]), "y": ("y", lats[:, 0])},
        )

        g.attrs = {
            "Coordinate System": cs,
            "alfori": alfori,
            "xori": xori,
            "yori": yori,
        }

        return g

    def to_file(self, filename, **kwargs):

        logger.info("writing grid to file {}".format(filename))

        with open(filename, "w") as f:
            f.write("Coordinate System= {}\n".format(self.Dataset.attrs["Coordinate System"]))
            f.write("{} {}\n".format(self.Dataset.lons.shape[1], self.Dataset.lons.shape[0]))
            f.write(
                "{} {} {}\n".format(
                    self.Dataset.attrs["xori"],
                    self.Dataset.attrs["yori"],
                    self.Dataset.attrs["alfori"],
                )
            )
            for i in range(self.Dataset.lons.shape[0]):
                f.write("ETA=  {} ".format(i + 1))
                f.write(" ".join(map(str, self.Dataset.lons[i, :].values)))
                f.write("\n")
            for i in range(self.Dataset.lats.shape[0]):
                f.write("ETA=  {} ".format(i + 1))
                f.write(" ".join(map(str, self.Dataset.lats[i, :].values)))
                f.write("\n")


class tri2d:
    """Unstructured triangular 2d mesh"""

    def __init__(self, **kwargs):

        mesh_file = kwargs.get("mesh_file", None)
        mesh_generator = kwargs.get("mesh_generator", None)
        geo = kwargs.get("geometry", None)
        coasts = kwargs.get("coastlines", None)
        boundary = kwargs.get("boundary", None)

        if geo == "global":
            kwargs.update({"gglobal": True})

        if mesh_file:

            self.Dataset = self.read_file(mesh_file)

        elif mesh_generator == "gmsh":

            if boundary is None:
                self.boundary = pb.get_boundaries(**kwargs)
            else:
                self.boundary = boundary

            g, bg = mgmsh.get(self.boundary.contours, **kwargs)  # create mesh with GMSH

            self.Dataset = g

            self.bgmesh = bg

        elif mesh_generator == "jigsaw":

            if boundary is None:
                self.boundary = pb.get_boundaries(**kwargs)
            else:
                self.boundary = boundary

            g, bg = mjigsaw.get(self.boundary.contours, **kwargs)  # create mesh with JIGSAW

            self.Dataset = g

            self.bgmesh = bg

        else:

            self.Dataset = None

    @staticmethod
    def read_file(hgrid, **kwargs):

        logger.info("read mesh file {}".format(hgrid))

        # read file
        df = pd.read_csv(hgrid, header=0, low_memory=False)
        df = df.dropna(axis=1)
        df.columns = ["data"]

        # extract number of elements, number of nodes
        ni, nj = df.iloc[0].str.split()[0]
        ni = int(ni)
        nj = int(nj)

        # read lon,lat,depth for all nodes
        q = pd.DataFrame(df.loc[1:nj, "data"].str.split().values.tolist())
        q = q.drop(q.columns[0], axis=1)
        q = q.apply(pd.to_numeric)
        #  q.reset_index(inplace=True, drop=True)
        q.columns = ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]
        q.index.name = "nSCHISM_hgrid_node"

        # create xarray of grid
        mesh = q.loc[:, ["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y"]].to_xarray()
        mesh = mesh.drop_vars("nSCHISM_hgrid_node")

        # create xarray of depth
        depth = q.loc[:, "depth"].to_xarray()
        depth = depth.drop_vars("nSCHISM_hgrid_node")

        # read connectivity
        e = pd.DataFrame(df.loc[nj + 1 : nj + ni, "data"].str.split().values.tolist())
        e = e.drop(e.columns[0], axis=1)
        e = e.apply(pd.to_numeric)

        ncolumns = e.loc[:, e.columns[0]].max()

        if ncolumns == 3:
            e.columns = ["nv", "a", "b", "c"]
        else:
            e.columns = ["nv", "a", "b", "c", "d"]

        e.loc[:, e.columns[1:]] = e.loc[:, e.columns[1:]].values - 1  # convert to python (index starts from 0)

        # create xarray of tessellation
        els = xr.DataArray(
            e.loc[:, e.columns[1:]].values,
            dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
            name="SCHISM_hgrid_face_nodes",
        )

        # Open boundaries
        n0 = df[df.data.str.contains("open boundaries")].index
        n0 = n0.values[0]
        nob = df.loc[n0, "data"].split("=")[0].strip()
        nob = int(nob)
        nobn = df.loc[n0 + 1, "data"].split("=")[0].strip()
        nobn = int(nobn)

        odic = {}
        ottr = []
        idx = n0 + 2
        for nl in range(nob):
            nn = df.loc[idx, "data"].split("=")[0].strip()
            nn = int(nn)
            label = df.loc[idx, "data"].split("=")[1]
            label = label[label.index("open") :]
            ottr.append([nn, label])
            nodes = df.loc[idx + 1 : idx + nn, "data"].astype("int")
            idx = idx + nn + 1

            oi = pd.DataFrame({"node": nodes, "type": "open", "id": int(label[-1])})
            odic.update({label: oi})

        try:
            dfo = pd.concat(odic).droplevel(0).reset_index(drop=True)
        except ValueError:
            dfo = pd.DataFrame()

        # Land boundaries
        n1 = df[df.data.str.contains("land boundaries")].index
        n1 = n1.values[0]

        nlb = df.loc[n1, "data"].split("=")[0].strip()
        nlb = int(nlb)

        nlbn = df.loc[n1 + 1, "data"].split("=")[0].strip()
        nlbn = int(nlbn)

        ldic = {}
        attr = []
        idx = n1 + 2
        ili = -1
        lli = 1001
        for nl in range(nlb):
            nn, etype = df.loc[idx, "data"].split("=")[0].strip().split(" ")
            nn = int(nn)
            etype = int(etype)
            label = df.loc[idx, "data"].split("=")[1]
            label = label[label.index("land") :]
            attr.append([nn, etype, label])
            nodes = df.loc[idx + 1 : idx + nn, "data"].astype(int)
            idx = idx + nn + 1

            li = pd.DataFrame({"node": nodes})
            tt = ["land" if etype == 0 else "island"]
            idi = [ili if etype == 1 else 1000 + int(label[-1])]
            li["type"] = tt[0]
            li["id"] = idi[0]
            ldic.update({label: li})
            if tt[0] == "land":
                lli += 1
            elif tt[0] == "island":
                ili -= 1
            else:
                raise ValueError(f"mesh boundaries error")
        try:
            dfl = pd.concat(ldic).droplevel(0).reset_index(drop=True)
        except ValueError:
            dfl = pd.DataFrame()

        # concat boundaries
        bbs = pd.concat([dfo, dfl])

        bbs.node = bbs.node - 1  # start_index = 0
        bbs = bbs.reset_index(drop=True)  # reset index
        bbs = bbs[["type", "node", "id"]]  # set column order
        bbs = bbs.sort_values(["type", "id", "node"]).reset_index(drop=True)  # sort
        bbs.index.name = "bnodes"

        # merge to one xarray DataSet
        g = xr.merge([mesh, depth, els, bbs.to_xarray()])

        g.attrs = {}

        return g

    def to_file(self, filename, **kwargs):

        logger.info("writing mesh to file {}".format(filename))

        nn = self.Dataset.SCHISM_hgrid_node_x.size
        n3e = self.Dataset.nSCHISM_hgrid_face.size

        with open(filename, "w") as f:
            f.write("\t uniform.gr3\n")
            f.write("\t {} {}\n".format(n3e, nn))

        q = self.Dataset[["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"]].to_dataframe()

        q.index = np.arange(1, len(q) + 1)

        q.to_csv(
            filename,
            index=True,
            sep="\t",
            header=None,
            mode="a",
            float_format="%.10f",
            columns=["SCHISM_hgrid_node_x", "SCHISM_hgrid_node_y", "depth"],
        )

        e = pd.DataFrame(
            self.Dataset.SCHISM_hgrid_face_nodes.dropna(dim="nMaxSCHISM_hgrid_face_nodes").values,
            columns=["a", "b", "c"],
        )

        e["nv"] = e.apply(lambda row: row.dropna().size, axis=1)

        e.index = np.arange(1, len(e) + 1)

        e = e.dropna(axis=1).astype(int)

        e.loc[:, ["a", "b", "c"]] = e.loc[:, ["a", "b", "c"]] + 1  # convert to fortran (index starts from 1)

        e.to_csv(
            filename,
            index=True,
            sep="\t",
            header=None,
            mode="a",
            columns=["nv", "a", "b", "c"],
        )

        bs = self.Dataset[["node", "id", "type"]].to_dataframe()

        # open boundaries
        number_of_open_boundaries = bs.loc[bs.type == "open"].id
        if not number_of_open_boundaries.empty:
            number_of_open_boundaries = number_of_open_boundaries.max()
        else:
            number_of_open_boundaries = 0
        number_of_open_boundaries_nodes = bs.loc[bs.type == "open"].shape[0]

        if number_of_open_boundaries > 0:
            with open(filename, "a") as f:

                f.write("{} = Number of open boundaries\n".format(number_of_open_boundaries))
                f.write("{} = Total number of open boundary nodes\n".format(number_of_open_boundaries_nodes))

                for i in range(1, number_of_open_boundaries + 1):
                    dat = bs.loc[bs.id == i, "node"] + 1  # fortran
                    f.write("{} = Number of nodes for open boundary {}\n".format(dat.size, i))
                    dat.to_csv(f, index=None, header=False)

        else:

            with open(filename, "a") as f:

                f.write("{} = Number of open boundaries\n".format(0))
                f.write("{} = Total number of open boundary nodes\n".format(0))

        # land boundaries

        number_of_land_boundaries = bs.loc[bs.type == "land"].id
        if not number_of_land_boundaries.empty:
            number_of_land_boundaries = number_of_land_boundaries.max() - 1000
        else:
            number_of_land_boundaries = 0
        number_of_land_boundaries_nodes = bs.loc[bs.type == "land"].shape[0]

        number_of_island_boundaries = bs.loc[bs.type == "island"].id
        if not number_of_island_boundaries.empty:
            number_of_island_boundaries = number_of_island_boundaries.min()
        else:
            number_of_island_boundaries = 0
        number_of_island_boundaries_nodes = bs.loc[bs.type == "island"].shape[0]

        nlb = number_of_land_boundaries - number_of_island_boundaries
        nlbn = number_of_land_boundaries_nodes + number_of_island_boundaries_nodes

        if nlb > 0:
            with open(filename, "a") as f:
                f.write("{} = Number of land boundaries\n".format(nlb))
                f.write("{} = Total number of land boundary nodes\n".format(nlbn))
                ik = 1
        else:
            with open(filename, "a") as f:
                f.write("{} = Number of land boundaries\n".format(0))
                f.write("{} = Total number of land boundary nodes\n".format(0))

        if number_of_land_boundaries > 0:
            with open(filename, "a") as f:

                for i in range(1001, 1000 + number_of_land_boundaries + 1):
                    dat_ = bs.loc[bs.id == i]
                    dat = dat_.node + 1  # fortran

                    f.write("{} {} = Number of nodes for land boundary {}\n".format(dat.size, 0, ik))
                    dat.to_csv(f, index=None, header=False)
                    ik += 1

        if number_of_island_boundaries < 0:

            with open(filename, "a") as f:

                for i in range(-1, number_of_island_boundaries - 1, -1):
                    dat_ = bs.loc[bs.id == i]
                    dat = dat_.node + 1  # fortran

                    f.write("{} {} = Number of nodes for land boundary {}\n".format(dat.size, 1, ik))
                    dat.to_csv(f, index=None, header=False)
                    ik += 1

    def validate(self, **kwargs):

        # ---------------------------------------------------------------------
        logger.info("start mesh validation\n")
        # ---------------------------------------------------------------------

        path = kwargs.get("rpath", "./")

        if not os.path.exists(path):
            os.makedirs(path)

        # save bctides.in
        bs = self.Dataset[["node", "id", "type"]].to_dataframe()
        # open boundaries
        number_of_open_boundaries = np.nan_to_num(bs.loc[bs.type == "open"].id.max()).astype(int)
        number_of_open_boundaries_nodes = bs.loc[bs.type == "open"].shape[0]

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

        # save params.nml
        config_file = DATA_PATH + "param_val.nml"
        params = f90nml.read(config_file)
        params.write(path + "param.nml", force=True)

        # save hgrid.gr3
        self.to_file(filename=path + "hgrid.gr3")

        try:
            bin_path = os.environ["SCHISM"]
        except:
            bin_path = kwargs.get("epath", None)

        if bin_path is None:
            # ------------------------------------------------------------------------------
            logger.warning("Schism executable path (epath) not given -> using default \n")
            # ------------------------------------------------------------------------------
            bin_path = "schism"

        tools.create_mpirun_script(
            target_dir=path,
            cmd=bin_path,
            script_name="launchSchism.sh",
            ncores=1,
        )

        # note that cwd is the folder where the executable is
        ex = subprocess.Popen(
            args=["./launchSchism.sh"],
            cwd=path,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )  # , bufsize=1)

        out, err = ex.communicate()[:]

        if "successfully" in str(out):

            # ---------------------------------------------------------------------
            logger.info("mesh is validated for SCHISM\n")
            # ---------------------------------------------------------------------
            return True
        else:
            logger.debug(str(out))
            # ---------------------------------------------------------------------
            logger.info("mesh fails.. exiting \n")
            # ---------------------------------------------------------------------
            return False

    def verify(self, **kwargs):

        shp = kwargs.get("coastlines", None)

        if shp is not None:
            r = verify(self, shp)
            return r
        else:
            logger.warning("No coastlines provided")
