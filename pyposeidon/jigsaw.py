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
from pyposeidon.utils.sort import *
import pyposeidon.dem as pdem
from pyposeidon.utils.hfun import *
from pyposeidon.utils.spline import *
from pyposeidon.utils.tag import *
import logging

logger = logging.getLogger("pyposeidon")


DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/misc/"
TEST_DATA_PATH = os.path.dirname(pyposeidon.__file__) + "/tests/data/"


def geo(df, path=".", tag="jigsaw"):

    fgeo = path + tag + "-geo.msh"
    # write header
    with open(fgeo, "w") as f:
        f.write("#{}; created by pyposeidon\n".format(tag + "-geo.msh"))
        f.write("MSHID=2;EUCLIDEAN-MESH\n")
        f.write("NDIMS=2\n")
        f.write("POINT={}\n".format(df.shape[0]))

    # write lines
    with open(fgeo, "a") as f:
        for line in df.index.levels[0]:
            df.loc[line].to_csv(f, index=False, header=0, columns=["lon", "lat", "z"], sep=";")

    edges = pd.DataFrame([])  # initiate

    # create edges
    for line in df.index.levels[0]:
        i0 = edges.shape[0]
        ie = df.loc[line].shape[0] + edges.shape[0]
        dline = df.loc[line].copy()
        dline.index = range(i0, ie)
        dline.loc[:, "ie"] = dline.index.values + 1
        dline.loc[dline.index[-1], "ie"] = i0
        dout = dline.reset_index().loc[:, ["index", "ie", "tag"]]

        dout["tag1"] = dline.loc[dline.ie.values, "tag"].values.astype(int)
        dout["que"] = np.where(((dout["tag"] != dout["tag1"]) & (dout["tag"] > 0)), dout["tag1"], dout["tag"])
        dout = dout.reset_index().loc[:, ["index", "ie", "que"]]
        edges = edges.append(dout)

    # write header
    with open(fgeo, "a") as f:
        f.write("EDGE2={}\n".format(edges.shape[0]))

    with open(fgeo, "a") as f:
        edges.to_csv(f, index=False, header=0, sep=";")


def read_msh(fmsh):

    grid = pd.read_csv(fmsh, header=0, names=["data"], index_col=None, low_memory=False)
    npoints = int(grid.loc[2].str.split("=")[0][1])

    nodes = pd.DataFrame(
        grid.loc[3 : 3 + npoints - 1, "data"].str.split(";").values.tolist(), columns=["lon", "lat", "z"]
    )

    ie = grid[grid.data.str.contains("EDGE")].index.values[0]
    nedges = int(grid.loc[ie].str.split("=")[0][1])
    edges = pd.DataFrame(
        grid.loc[ie + 1 : ie + nedges, "data"].str.split(";").values.tolist(), columns=["e1", "e2", "e3"]
    )

    i3 = grid[grid.data.str.contains("TRIA")].index.values[0]
    ntria = int(grid.loc[i3].str.split("=")[0][1])
    tria = pd.DataFrame(
        grid.loc[i3 + 1 : i3 + ntria + 1, "data"].str.split(";").values.tolist(), columns=["a", "b", "c", "d"]
    )

    return [nodes, edges, tria]


def jigsaw(**kwargs):

    logger.info("Creating grid with JIGSAW\n")

    geometry = kwargs.get("geometry", None)

    # Set bgmesh
    bgmesh = kwargs.get("bgmesh", None)

    if bgmesh == "auto":

        try:

            logger.info("Read DEM")
            dem = pdem.dem(**kwargs)

            res_min = kwargs.get("resolution_min", 0.01)
            res_max = kwargs.get("resolution_max", 0.5)
            dhdx = kwargs.get("dhdx", 0.15)

            logger.info("Evaluate bgmesh")
            w = hfun(
                dem.Dataset.elevation, resolution_min=res_min, resolution_max=res_max, dhdx=dhdx
            )  # resolution in lat/lon degrees

            path = kwargs.get("rpath", "./bgmesh/")

            if not os.path.exists(path):  # check if run folder exists
                os.makedirs(path)

            logger.info("Save bgmesh to {}bgmesh.nc".format(path))
            w.to_netcdf(path + "bgmesh.nc")  # save bgmesh

            kwargs.update({"bgmesh": path + "bgmesh.nc"})

        except:

            logger.warning("bgmesh failed... continuing without background mesh size")

    if isinstance(geometry, dict):

        df, bmindx = tag_(**kwargs)

        gr = jigsaw_(df, bmindx, **kwargs)

    elif isinstance(geometry, str):

        if geometry == "global":

            bgmesh = hfun_(kwargs.get("coastlines", None), kwargs.get("res", 0.1), kwargs.get("R", 1.0))

            kwargs.update({"bgmesh": bgmesh})

            df = sgl(**kwargs)

            bmindx = df.tag.min()

            gr = jigsaw_(df, bmindx, **kwargs)

            convert = kwargs.get("to_lat_lon", True)

            if convert:
                # convert to lat/lon
                u, v = gr.SCHISM_hgrid_node_x.values, gr.SCHISM_hgrid_node_y.values

                rlon, rlat = to_lat_lon(u, v, R=kwargs.get("R", 1.0))

                gr["SCHISM_hgrid_node_x"].values = rlon
                gr["SCHISM_hgrid_node_y"].values = rlat

        else:

            df = jcustom(**kwargs)

            bmindx = df.tag.min()

            gr = jigsaw_(df, bmindx, **kwargs)

    return gr


def sgl(**kwargs):

    try:
        geo = gp.GeoDataFrame.from_file(kwargs.get("coastlines", None))
    except:
        logger.warning("coastlines is not a file, trying with geopandas Dataset")
        try:
            geo = kwargs.get("coastlines", None)
        except:
            logger.error("coastlines argument not valid ")
            sys.exit(1)

        # Manage coastlines
        logger.info("preparing coastlines")

        # ANTARTICA
        anta_mask = geo.bounds.miny < geo.bounds.miny.min() + 0.1  # indentify antartica
        anta = geo.loc[anta_mask]
        indx = anta.index  # keep index

        try:
            anta = pd.DataFrame(anta.boundary.values[0].coords[:], columns=["lon", "lat"])
        except:
            anta = pd.DataFrame(
                anta.boundary.explode()[0].coords[:], columns=["lon", "lat"]
            )  # convert boundary values to pandas
        d1 = anta.where(anta.lon == anta.lon.max()).dropna().index[1:]  # get artificial boundaries as -180/180
        d2 = anta.where(anta.lon == anta.lon.min()).dropna().index[1:]
        anta = anta.drop(d1).drop(d2)  # drop the points
        d3 = anta.where(anta.lat == anta.lat.min()).dropna().index  # drop lat=-90 line
        anta = anta.drop(d3)
        an = gp.GeoDataFrame(
            {
                "geometry": [shapely.geometry.LineString(anta.values)],
                "length": shapely.geometry.LineString(anta.values).length,
            },
            index=indx,
        )  # put together a LineString
        geo.loc[indx] = shapely.geometry.LineString(anta.values)  # put it back to geo

        # International Meridian
        m1 = geo[geo.bounds.minx == -180.0].index
        m2 = geo[geo.bounds.maxx == 180.0].index
        mm = np.concatenate((m1, m2))  # join them
        mm = [j for j in mm if j != indx]  # subtract antartica

        # convert to u,v (stereographic coordinates)
        for idx, poly in geo.iterrows():
            geo.loc[idx, "geometry"] = shapely.ops.transform(
                lambda x, y, z=None: to_stereo(x, y, R=kwargs.get("R", 1.0)), poly.geometry
            )

        w = geo.drop(indx)  # get all polygons
        ww = w.loc[mm]  # join the split polygons
        gw = gp.GeoDataFrame(
            geometry=list(ww.buffer(0.0001).unary_union)
        )  # merge the polygons that are split (around -180/180)

        w = w.drop(mm)
        # Check antartica LineString
        if not geo.iloc[indx].geometry.values[0].is_ring:
            ca = gp.GeoDataFrame(geometry=[shapely.geometry.LinearRing(geo.loc[indx].geometry.values[0])], index=indx)
            ca["geometry"] = shapely.geometry.LineString(ca.geometry.values[0])
        else:
            ca = geo.loc[indx]

        # PUT ALL TOGETHER
        geo = pd.concat([w, gw, ca], ignore_index=True).reset_index(drop=True)

    logger.info("storing boundaries")

    geo["tag"] = -(geo.index + 1)

    idx = 0
    dic = {}
    for i, line in geo.iloc[:-1].iterrows():
        lon = []
        lat = []
        try:
            for x, y in line.geometry.boundary.coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({"line{}".format(idx): {"lon": lon, "lat": lat, "tag": line.tag}})
            idx += 1
        except:
            for x, y in line.geometry.boundary[0].coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({"line{}".format(idx): {"lon": lon, "lat": lat, "tag": line.tag}})
            idx += 1

    for i, line in geo.iloc[-1:].iterrows():
        lon = []
        lat = []
        for x, y in line.geometry.coords[:]:
            lon.append(x)
            lat.append(y)
            dic.update({"line{}".format(idx): {"lon": lon, "lat": lat, "tag": line.tag}})

    dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}

    df = pd.concat(dict_of_df, axis=0)

    df["z"] = 0
    df = df.drop_duplicates()  # drop the repeat value on closed boundaries

    return df


def jcustom(**kwargs):

    geometry = kwargs.get("geometry", None)

    try:
        geo = gp.GeoDataFrame.from_file(geometry)
    except:
        logger.error("geometry argument not a valid file")
        sys.exit(1)

    idx = 0
    dic = {}
    for idx, line in geo.iterrows():
        lon = []
        lat = []
        for x, y in line.geometry.coords[:]:
            lon.append(x)
            lat.append(y)
        dic.update({"line{}".format(idx): {"lon": lon, "lat": lat, "tag": line.tag}})
        idx += 1

    dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}

    df = pd.concat(dict_of_df, axis=0)

    df["z"] = 0
    df = df.drop_duplicates()  # drop the repeat value on closed boundaries

    out_b = []
    for line in df.index.levels[0]:
        out_b.append(shapely.geometry.LineString(df.loc[line, ["lon", "lat"]].values))

    merged = shapely.ops.linemerge(out_b)
    merged = pd.DataFrame(merged.coords[:], columns=["lon", "lat"])
    merged = merged.drop_duplicates()
    match = df.drop_duplicates(["lon", "lat"]).droplevel(0)
    match = match.reset_index(drop=True)

    df1 = merged.sort_values(["lon", "lat"])
    df2 = match.sort_values(["lon", "lat"])
    df2.index = df1.index
    final = df2.sort_index()
    final = pd.concat([final], keys=["line0"])

    return final


def jigsaw_(df, bmindx, **kwargs):

    logger.info("Creating JIGSAW files\n")

    tag = kwargs.get("tag", "jigsaw")
    rpath = kwargs.get("rpath", ".")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    path = rpath + "/jigsaw/"
    if not os.path.exists(path):
        os.makedirs(path)

    geo(df, path=path, tag=tag)

    bgmesh = kwargs.get("bgmesh", None)

    if bgmesh is not None:

        if isinstance(bgmesh, str):
            dh = xr.open_dataset(bgmesh)
            to_hfun_grid(dh, path + tag + "-hfun.msh")  # write bgmesh file
        else:
            to_hfun_mesh(bgmesh, path + tag + "-hfun.msh")

    # write jig file
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

    execute = kwargs.get("execute_jigsaw", True)

    if execute:

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

        gr = to_dataset(bmindx, **kwargs)

        logger.info("..done creating mesh\n")

        return gr


def to_dataset(bmindx, **kwargs):

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

    # LAND Boundaries (negative tag)
    ib = 1
    ib_ = 1
    isl = []
    for ik in range(edges.e3.min(), 0):
        bb = np.unique(edges.loc[edges.e3 == ik, ["e1", "e2"]].values.flatten())
        bf = pd.concat([nodes.loc[bb]])
        bf["idx"] = bb
        if ik >= bmindx:
            flag = 0
            bf["indx"] = "land_boundary_{}_{}".format(ib, flag)
            c1 = True
            c2 = False
        else:
            flag = 1
            bf["indx"] = "land_boundary_{}_{}".format(ib_, flag)
            c1 = False
            c2 = True

        bf["flag"] = flag
        bf = bf.reset_index().set_index(["indx", "index"])
        bf.index.names = ["", ""]

        if not bf.empty:
            isl.append(bf)
            if c1 is True:
                ib += 1
            if c2 is True:
                ib_ += 1

    if isl:

        landb = pd.concat(isl)

        land_i = sorted(landb.index.levels[0], key=lambda x: [int(x.split("_")[3]), int(x.split("_")[2])])

    else:

        land_i = []

    # WATER Boundaries (positive tag)

    wb = 1
    wbs = []
    for ik in range(1, edges.e3.max() + 1):
        bb = np.unique(edges.loc[edges.e3 == ik, ["e1", "e2"]].values.flatten())
        bf = pd.concat([nodes.loc[bb]], keys=["open_boundary_{}".format(wb)])
        bf["idx"] = bb

        wbs.append(bf)
        wb += 1

    if wbs:

        openb = pd.concat(wbs)

        open_i = sorted(openb.index.levels[0], key=lambda x: int(x.split("_")[2]))

        if openb.index.levels[0].shape[0] == 1:  # sort the nodes if open box

            pts = openb[["lon", "lat"]].values

            origin = [openb.mean()["lon"], openb.mean()["lat"]]

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

        open_i = []

    # Fix an issue with jumping identifiers

    nns = np.array([int(x.split("_")[-1]) for x in open_i])
    try:
        if nns.max() != nns.size:
            nns = np.arange(nns.size) + 1
            open_i = ["open_boundary_{}".format(x) for x in nns]
            lidx = dict(zip(openb.index.levels[0].to_list(), open_i))
            openb.rename(index=lidx, level=0, inplace=True)
    except:
        pass  # continue is nns=[]

    # MAKE Dataset

    els = xr.DataArray(
        tria.loc[:, ["a", "b", "c"]].values,
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        nodes.loc[:, ["lon", "lat"]]
        .to_xarray()
        .rename({"index": "nSCHISM_hgrid_node", "lon": "SCHISM_hgrid_node_x", "lat": "SCHISM_hgrid_node_y"})
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    dep = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    # open boundaries
    op = []
    nop = []
    o_label = []
    for line in open_i:
        op.append(openb.loc[line, "idx"].values)
        nop.append(openb.loc[line, "idx"].size)
        o_label.append(line)

    oattr = pd.DataFrame({"label": o_label, "nps": nop})
    oattr["type"] = np.nan
    oattr.set_index("label", inplace=True, drop=True)

    obns = [j for i in op for j in i]
    odf = pd.DataFrame({"node": obns, "type": 0, "id": 0}, index=np.arange(len(obns)))
    idx = 1
    for l in range(len(op)):
        odf.loc[odf.node.isin(op[l]), "id"] = idx
        odf.loc[odf.node.isin(op[l]), "type"] = oattr.iloc[l].type
        idx += 1

    # land boundaries

    lp = []
    nlp = []
    l_label = []
    itype = []
    for line in land_i:
        lp.append(landb.loc[line, "idx"].values)
        nlp.append(landb.loc[line, "idx"].size)
        itype.append(landb.loc[line, "flag"].values[0])
        l_label.append(line)

    lattr = pd.DataFrame({"label": l_label, "nps": nlp})
    lattr["type"] = itype
    lattr.set_index("label", inplace=True, drop=True)

    lbns = [j for i in lp for j in i]
    ldf = pd.DataFrame({"node": lbns, "type": 1, "id": 0}, index=np.arange(len(lbns)))

    idx = -1
    for l in range(len(lp)):
        ldf.loc[ldf.node.isin(lp[l]), "id"] = idx
        ldf.loc[ldf.node.isin(lp[l]), "type"] = lattr.iloc[l].type
        idx -= 1

    if not ldf.empty:
        lk = -1
        for k in range(-1, ldf.id.min() - 1, -1):
            if not ldf.loc[ldf.id == k].empty:
                ldf.loc[ldf.id == k, "id"] = lk
                lk -= 1

    tbf = pd.concat([odf, ldf])
    tbf.index.name = "bnodes"

    gr = xr.merge([nod, dep, els, tbf.to_xarray()])  # total

    return gr
