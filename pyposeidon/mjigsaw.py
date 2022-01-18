"""
Jigsaw module

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import pandas as pd
import numpy as np
import geopandas as gp
import xarray as xr
import os
import shapely
import subprocess
from tqdm import tqdm
import sys

from pyposeidon.utils.stereo import to_lat_lon, to_stereo
from pyposeidon.utils.global_bgmesh import make_bgmesh_global
from pyposeidon.utils.sort import *
import pyposeidon.dem as pdem
from pyposeidon.utils.hfun import *
from pyposeidon.utils.spline import *
from pyposeidon.utils.tag import *
import logging

logger = logging.getLogger(__name__)


DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/misc/"
TEST_DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/tests/data/"


def to_geo(df, path=".", tag="jigsaw"):

    fgeo = path + tag + "-geo.msh"
    # write header
    with open(fgeo, "w") as f:
        f.write("#{}; created by pyposeidon\n".format(tag + "-geo.msh"))
        f.write("MSHID=2;EUCLIDEAN-MESH\n")
        f.write("NDIMS=2\n")
        f.write("POINT={}\n".format(df.nps.sum()))

    # outer contour
    df_ = df.loc[df.tag != "island"].reset_index(drop=True)  # all external contours

    if not df_.empty:

        # store xy in a DataFrame
        dic = {}
        for k, d in df_.iterrows():
            out = pd.DataFrame(d.geometry.coords[:], columns=["x", "y"])
            out["lindex"] = d.lindex
            out = out.drop_duplicates(["x", "y"])
            if d.tag == "land":  # drop end points in favor or open tag
                out = out[1:-1]
            dic.update({"line{}".format(k): out})
        o1 = pd.concat(dic, axis=0).droplevel(0).reset_index(drop=True)
        o1 = o1.drop_duplicates(["x", "y"])

        # Do linemerge of outer contours
        lss = df_.geometry.values
        merged = shapely.ops.linemerge(lss)
        o2 = pd.DataFrame({"x": merged.xy[0], "y": merged.xy[1]})  # convert to DataFrame
        o2 = o2.drop_duplicates()

        outer = o2.merge(o1)  # merge to transfer the lindex

        # write lines & compute edges
        # outer to file
        with open(fgeo, "a") as f:
            outer["z"] = 0
            outer = outer.drop_duplicates(["x", "y"])
            # nodes
            outer.to_csv(f, index=False, header=0, columns=["x", "y", "z"], sep=";")
            # compute edges

        edges = [
            list(a)
            for a in zip(
                np.arange(outer.shape[0]),
                np.arange(outer.shape[0]) + 1,
                outer.lindex.values,
            )
        ]  # outer
        edges[-1][1] = 0

        # sort (duplicated bounds)
        edges = pd.DataFrame(edges, columns=["index", "ie", "lindex"])

        edges["lindex1"] = (
            pd.concat([outer.loc[1:, "lindex"], outer.loc[0:0, "lindex"]]).reset_index(drop=True).astype(int)
        )
        edges["que"] = np.where(
            ((edges["lindex"] != edges["lindex1"]) & (edges["lindex"] > 0) & (edges["lindex"] < 1000)),
            edges["lindex1"],
            edges["lindex"],
        )

        edges = edges.reset_index().loc[:, ["index", "ie", "que"]]

        edges = edges.values.tolist()

    else:

        edges = []

    # the rest to file
    with open(fgeo, "a") as f:
        for k, d in df.loc[df.tag == "island"].iterrows():
            out = pd.DataFrame(d.geometry.coords[:], columns=["x", "y"])
            out["z"] = 0
            out = out.drop_duplicates(["x", "y"])

            # nodes
            out.to_csv(f, index=False, header=0, columns=["x", "y", "z"], sep=";")
            # compute edges
            i0 = len(edges)
            ie = out.shape[0] + len(edges)

            lindex = d.lindex

            for m in range(i0, ie):
                edges.append([m, m + 1, lindex])

            edges[-1][1] = i0

    edges = pd.DataFrame(edges)  # convert to pandas

    # write header
    with open(fgeo, "a") as f:
        f.write("EDGE2={}\n".format(edges.shape[0]))

    with open(fgeo, "a") as f:
        edges.to_csv(f, index=False, header=0, sep=";")


def read_msh(fmsh):

    grid = pd.read_csv(fmsh, header=0, names=["data"], index_col=None, low_memory=False)
    npoints = int(grid.loc[2].str.split("=")[0][1])

    nodes = pd.DataFrame(
        grid.loc[3 : 3 + npoints - 1, "data"].str.split(";").values.tolist(),
        columns=["x", "y", "z"],
    )

    ie = grid[grid.data.str.contains("EDGE")].index.values[0]
    nedges = int(grid.loc[ie].str.split("=")[0][1])
    edges = pd.DataFrame(
        grid.loc[ie + 1 : ie + nedges, "data"].str.split(";").values.tolist(),
        columns=["e1", "e2", "e3"],
    )

    i3 = grid[grid.data.str.contains("TRIA")].index.values[0]
    ntria = int(grid.loc[i3].str.split("=")[0][1])
    tria = pd.DataFrame(
        grid.loc[i3 + 1 : i3 + ntria + 1, "data"].str.split(";").values.tolist(),
        columns=["a", "b", "c", "d"],
    )

    return [nodes, edges, tria]


def make_bgmesh(contours, **kwargs):

    gglobal = kwargs.get("gglobal", False)

    if gglobal:
        lon_min = -180.0
        lon_max = 180.0
        lat_min = -90.0
        lat_max = 90.0
    else:
        lon_min = contours.bounds.minx.min()
        lon_max = contours.bounds.maxx.max()
        lat_min = contours.bounds.miny.min()
        lat_max = contours.bounds.maxy.max()

    kwargs_ = kwargs.copy()
    kwargs_.pop("lon_min", None)
    kwargs_.pop("lon_max", None)
    kwargs_.pop("lat_min", None)
    kwargs_.pop("lat_max", None)

    dem = kwargs.get("dem_source", None)

    if not isinstance(dem, xr.Dataset):
        logger.info("Read DEM")
        dem = pdem.Dem(lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, **kwargs_)
        dem = dem.Dataset

    res_min = kwargs_.pop("resolution_min", 0.01)
    res_max = kwargs_.pop("resolution_max", 0.5)
    dhdx = kwargs_.pop("dhdx", 0.15)

    rpath = kwargs.get("rpath", ".")
    tag = kwargs.get("tag", "jigsaw")
    fpos = rpath + "/jigsaw/" + tag + "-hfun.msh"

    if gglobal:

        logger.info("Evaluate global bgmesh")
        nds, lms = make_bgmesh_global(contours, fpos, dem, **kwargs)
        logger.info("Saving global background scale file")
        dh = to_global_hfun(nds, lms, fpos, **kwargs)

    else:

        logger.info("Evaluate bgmesh")
        dh = get_hfun(
            dem, resolution_min=res_min, resolution_max=res_max, dhdx=dhdx, **kwargs_
        )  # resolution in lat/lon degrees

        logger.info("Saving background scale file")
        to_hfun_grid(dh, fpos)  # write bgmesh file

    return dh


def to_dataset(**kwargs):

    logger.info("..reading mesh\n")

    tag = kwargs.get("tag", "jigsaw")
    rpath = kwargs.get("rpath", ".")

    [nodes, edges, tria] = read_msh(rpath + "/jigsaw/" + tag + ".msh")
    nodes = nodes.apply(pd.to_numeric)
    tria = tria.apply(pd.to_numeric)
    edges = edges.apply(pd.to_numeric)

    # Look for hanging nodes
    tri3 = tria.values[:, :3]
    q = np.unique(tri3.flatten())  # all the unique nodes in elements

    dq = list(set(range(nodes.shape[0])) - set(q))  # the ones that are in gcrop but not in elems

    dq.sort()
    nodes = nodes.drop(dq)  # drop nodes
    nodes = nodes.rename_axis("tag").reset_index()  # reset index

    ### Re-index tessalation

    A, idxA = np.unique(nodes["tag"], return_inverse=True)
    B, idxB = np.unique(tria["a"], return_inverse=True)
    IDX = np.in1d(A, B)
    tria["a"] = idxA[IDX][idxB]
    B, idxB = np.unique(tria["b"], return_inverse=True)
    IDX = np.in1d(A, B)
    tria["b"] = idxA[IDX][idxB]
    B, idxB = np.unique(tria["c"], return_inverse=True)
    IDX = np.in1d(A, B)
    tria["c"] = idxA[IDX][idxB]

    # Drop invalid edges
    drop_e = edges.loc[edges.e1.isin(dq) | edges.e2.isin(dq)].index
    edges = edges.drop(drop_e).reset_index(drop=True)
    ### Re-index edges
    A, idxA = np.unique(nodes["tag"], return_inverse=True)
    B, idxB = np.unique(edges["e1"], return_inverse=True)
    IDX = np.in1d(A, B)
    edges["e1"] = idxA[IDX][idxB]
    B, idxB = np.unique(edges["e2"], return_inverse=True)
    IDX = np.in1d(A, B)
    edges["e2"] = idxA[IDX][idxB]
    # clean up
    nodes = nodes.drop("tag", axis=1)

    # Boundaries

    # LAND Boundaries (1000+ tag)
    lnd = []
    for ik in range(1001, edges.e3.max() + 1):
        bb = np.unique(edges.loc[edges.e3 == ik, ["e1", "e2"]].values.flatten())
        bf = pd.concat([nodes.loc[bb]], keys=["land_boundary_{}".format(ik - 1000)])
        bf["id"] = ik

        if not bf.empty:
            lnd.append(bf)

    if lnd:
        landb = pd.concat(lnd)
    else:
        landb = pd.DataFrame([])

    # ISLAND Boundaries (negative tag)
    isl = []

    for ik in range(edges.e3.min(), 0):
        bb = np.unique(edges.loc[edges.e3 == ik, ["e1", "e2"]].values.flatten())
        bf = pd.concat([nodes.loc[bb]], keys=["island_boundary_{}".format(-ik)])
        bf["id"] = ik

        if not bf.empty:
            isl.append(bf)

    if isl:
        islandb = pd.concat(isl)
    else:
        islandb = pd.DataFrame([])

    # Open Boundaries (positive tag < 1000)
    nr_open = edges.loc[(edges.e3 > 0) & (edges.e3 < 1000)].e3
    if not nr_open.empty:
        nr_open = nr_open.max()
    else:
        nr_open = 0

    wbs = []
    for ik in range(0, nr_open + 1):
        bb = np.unique(edges.loc[edges.e3 == ik, ["e1", "e2"]].values.flatten())
        bf = pd.concat([nodes.loc[bb]], keys=["open_boundary_{}".format(ik)])
        bf["id"] = ik

        wbs.append(bf)

    if wbs:

        openb = pd.concat(wbs)

        if openb.index.levels[0].shape[0] == 1:  # sort the nodes if open box

            pts = openb[["x", "y"]].values

            origin = [openb.mean()["x"], openb.mean()["y"]]

            refvec = [1, 0]

            sps = sorted(pts, key=lambda po: clockwiseangle_and_distance(po, origin, refvec))

            sps = np.array(sps)

            # reorder openb
            fs = [tuple(lst) for lst in sps]
            ss = [tuple(lst) for lst in pts]

            index_dict = dict((value, idx) for idx, value in enumerate(ss))
            idx = [index_dict[x] for x in fs]

            openb = openb.iloc[idx]

    else:

        openb = pd.DataFrame([])

    # convert if global
    gglobal = kwargs.get("gglobal", False)
    if gglobal:
        R = kwargs.get("R", 1.0)
        rlon, rlat = to_lat_lon(nodes.x, nodes.y, R=R)
        nodes["x"] = rlon
        nodes["y"] = rlat

    # MAKE Dataset

    els = xr.DataArray(
        tria.loc[:, ["a", "b", "c"]].values,
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        nodes.loc[:, ["x", "y"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "x": "SCHISM_hgrid_node_x",
                "y": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    dep = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    # open boundaries
    if not openb.empty:
        odf = openb.reset_index()[["level_0", "level_1", "id"]]
        odf.columns = ["type", "node", "id"]
        odf.type = [x.split("_")[0] for x in odf.type]
    else:
        odf = None

    # land boundaries
    if not landb.empty:
        ldf = landb.reset_index()[["level_0", "level_1", "id"]]
        ldf.columns = ["type", "node", "id"]
        ldf.type = [x.split("_")[0] for x in ldf.type]
    else:
        ldf = None

    # island boundaries
    if not islandb.empty:
        idf = islandb.reset_index()[["level_0", "level_1", "id"]]
        idf.columns = ["type", "node", "id"]
        idf.type = [x.split("_")[0] for x in idf.type]
    else:
        idf = None

    tbf = pd.concat([odf, ldf, idf])
    tbf = tbf.reset_index(drop=True)
    tbf.index.name = "bnodes"

    gr = xr.merge([nod, dep, els, tbf.to_xarray()])  # total

    return gr


def get(contours, **kwargs):

    logger.info("Creating JIGSAW files\n")

    tag = kwargs.get("tag", "jigsaw")
    rpath = kwargs.get("rpath", ".")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    path = rpath + "/jigsaw/"
    if not os.path.exists(path):
        os.makedirs(path)

    # HFUN FILE
    bgmesh = kwargs.get("bgmesh", None)

    if bgmesh is None:
        dem_source = kwargs.get("dem_source", None)
        if dem_source:
            bgmesh = "auto"
            kwargs.update({"bgmesh": "auto"})

    if bgmesh is not None:

        logger.info("Set background scale")

        if bgmesh.endswith(".nc"):

            try:
                dh = xr.open_dataset(bgmesh)

                if "longitude" in dh.coords:
                    to_hfun_grid(dh, path + tag + "-hfun.msh")  # write bgmesh file
                else:
                    to_hfun_mesh(dh, path + tag + "-hfun.msh")
            except:
                logger.warning("bgmesh failed... continuing without background mesh size")
                bgmesh = None

        elif bgmesh == "auto":

            dh = make_bgmesh(contours, **kwargs)

        elif bgmesh.endswith(".msh"):

            pass

    try:
        bg = dh
    except:
        bg = None

    # GEO FILE
    to_geo(contours, path=path, tag=tag)

    # JIG FILE
    fjig = path + "/" + tag + ".jig"

    with open(fjig, "w") as f:
        f.write("GEOM_FILE ={}\n".format(tag + "-geo.msh"))
        f.write("MESH_FILE ={}\n".format(tag + ".msh"))
        if bgmesh:
            f.write("HFUN_FILE ={}\n".format(tag + "-hfun.msh"))
        f.write("HFUN_SCAL = ABSOLUTE\n")
        f.write("HFUN_HMAX = Inf\n")
        f.write("HFUN_HMIN = 0.0\n")
        f.write("MESH_DIMS = 2\n")
        f.write("MESH_TOP1 = TRUE\n")
        #        f.write('MESH_TOP2 = TRUE\n')
        f.write("MESH_EPS1 = 1.0\n")
        f.write("MESH_RAD2 = 1\n")
        f.write("GEOM_FEAT = TRUE\n")
        f.write("VERBOSITY = 2")

    calc_dir = rpath + "/jigsaw/"

    # EXECUTE
    setup_only = kwargs.get("setup_only", False)

    if not setup_only:

        # ---------------------------------
        logger.info("executing jigsaw\n")
        # ---------------------------------

        # execute jigsaw
        ex = subprocess.Popen(
            args=["jigsaw {}".format(tag + ".jig")],
            cwd=calc_dir,
            shell=True,
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )  # , bufsize=1)

        with open(calc_dir + "err.log", "w") as f:
            for line in iter(ex.stderr.readline, b""):
                f.write(line.decode(sys.stdout.encoding))
        #        logger.info(line.decode(sys.stdout.encoding))
        ex.stderr.close()

        with open(calc_dir + "run.log", "w") as f:
            for line in iter(ex.stdout.readline, b""):
                f.write(line.decode(sys.stdout.encoding))
        #        logger.info(line.decode(sys.stdout.encoding))
        ex.stdout.close()

        # ---------------------------------
        logger.info("Jigsaw FINISHED\n")
        # ---------------------------------

        gr = to_dataset(**kwargs)

        logger.info("..done creating mesh\n")

        return gr, bg

    else:

        gr = None

        return gr, bg
