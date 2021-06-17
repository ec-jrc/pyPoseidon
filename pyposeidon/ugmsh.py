"""
gmsh module

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
from tqdm import tqdm
import sys
import gmsh

from pyposeidon.utils.spline import *
from pyposeidon.utils.stereo import to_lat_lon
from pyposeidon.utils.pos import *
import pyposeidon.dem as pdem

import logging

logger = logging.getLogger("pyposeidon")


def read_gmsh(mesh, **kwargs):

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()

    gmsh.open(mesh)

    logger.info("Analyze grid")

    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    nodes = pd.DataFrame(coord.reshape(-1, 3), columns=["x", "y", "z"])

    elementTags2, nodeTags2 = gmsh.model.mesh.getElementsByType(2)

    elems = nodeTags2.reshape(-1, 3)

    tria = pd.DataFrame(elems - 1, columns=["a", "b", "c"])

    # boundaries

    bounds = []

    bgs = gmsh.model.getPhysicalGroups()

    bgs = pd.DataFrame(bgs[:-1], columns=["dim", "tag"])

    # open boundaries
    logger.info("open boundaries")

    obs = bgs.loc[bgs.tag < 1000]

    for row in obs.itertuples(index=True, name="Pandas"):
        onodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))

        db = pd.DataFrame({"node": onodes - 1})
        db["type"] = "open"
        db["id"] = getattr(row, "Index") + 1

        bounds.append(db)

    # land boundaries type
    logger.info("land boundaries")

    lbs = bgs.loc[(bgs.tag > 1000) & (bgs.tag < 2000)]
    lbs.reset_index(inplace=True, drop=True)

    for row in lbs.itertuples(index=True, name="Pandas"):
        lnodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))

        db = pd.DataFrame({"node": lnodes - 1})
        db["type"] = "land"
        db["id"] = 1000 + (getattr(row, "Index") + 1)

        bounds.append(db)

    # islands
    logger.info("islands")

    ibs = bgs.loc[bgs.tag > 2000]
    ibs.reset_index(inplace=True, drop=True)

    for row in ibs.itertuples(index=True, name="Pandas"):

        inodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))
        db = pd.DataFrame({"node": inodes - 1})
        db["type"] = "island"
        db["id"] = -(getattr(row, "Index") + 1)

        bounds.append(db)

    if bounds != []:
        bnodes = pd.concat(bounds).reset_index(drop=True)

        bnodes.index.name = "bnodes"

        bnodes = bnodes.drop_duplicates("node")

        bnodes["id"] = bnodes.id.astype(int)

    else:
        bnodes = pd.DataFrame({})

    # check if global and reproject
    sproj = kwargs.get("gglobal", False)
    if sproj:  # convert to lat/lon
        xd, yd = to_lat_lon(nodes.x, nodes.y)
        nodes["x"] = xd
        nodes["y"] = yd

    grid = pd.DataFrame({"lon": nodes.x, "lat": nodes.y})

    tri3 = tria.values

    logger.info("Finalize Dataset")

    ## make dataset
    els = xr.DataArray(
        tri3, dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"], name="SCHISM_hgrid_face_nodes"
    )

    nod = (
        grid.loc[:, ["lon", "lat"]]
        .to_xarray()
        .rename({"index": "nSCHISM_hgrid_node", "lon": "SCHISM_hgrid_node_x", "lat": "SCHISM_hgrid_node_y"})
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    dep = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    gr = xr.merge([nod, els, dep, bnodes.to_xarray()])

    gmsh.finalize()

    return gr


def gmsh_(df, **kwargs):

    logger.info("Creating grid with GMSH\n")

    rpath = kwargs.get("rpath", ".")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    gpath = os.path.join(rpath, "gmsh")
    if not os.path.exists(gpath):
        os.makedirs(gpath)

    make_gmsh(df, **kwargs)

    gr = read_gmsh(rpath + "/gmsh/mymesh.msh")

    return gr


def gset(df, **kwargs):

    logger.info("interpolate coastal points")

    lc = kwargs.get("lc", 0.5)

    df["lc"] = lc
    df = df.apply(pd.to_numeric)

    # Resample to equidistant points

    conts = np.unique(df.index[df.tag < 0].get_level_values(0))
    conts = [x for x in conts if x not in ["line0"]]  # except the outer LineString

    ibs = len(conts)

    ndfsfs = {}
    for ic in tqdm(range(ibs)):
        contour = conts[ic]
        curve = df.loc[contour, ["lon", "lat"]]
        curve = pd.concat([curve, curve.loc[0:0]]).reset_index(drop=True)
        di = spline(curve, ds=0.01, method="slinear")
        di["z"] = df.loc[contour].z.values[0]
        di["tag"] = df.loc[contour].tag.values[0].astype(int)
        di["lc"] = df.loc[contour].lc.values[0]
        ndfsfs.update({contour: di.drop_duplicates(["lon", "lat"])})

    df_ = pd.concat(ndfsfs, axis=0)
    df_["z"] = df_.z.values.astype(int)
    df_["tag"] = df_.tag.values.astype(int)

    # Line0

    logger.info("set outermost boundary")

    df0 = df.loc["line0"]

    mtag = df0.tag.min()
    mtag = mtag.astype(int)

    nd0 = {}
    for ic in tqdm(range(mtag, 0)):
        contour = df0.tag == ic
        curve = df0.loc[contour, ["lon", "lat"]].reset_index(drop=True)
        #    curve = pd.concat([curve,curve.loc[0:0]]).reset_index(drop=True)
        di = spline(curve, ds=0.01, method="slinear")
        di["z"] = df0.loc[contour].z.values[0]
        di["tag"] = df0.loc[contour].tag.values[0].astype(int)
        di["lc"] = df0.loc[contour].lc.values[0]
        nd0.update({ic: di.drop_duplicates(["lon", "lat"])})

    # Join Line0
    df0_ = df0.copy()
    for l in range(mtag, 0):
        #    print(l)
        idx = df0_.loc[df0_.tag == l].index
        df0_ = pd.concat([df0_.iloc[: idx[0]], nd0[l], df0_.iloc[idx[-1] + 1 :]])
        df0_.reset_index(drop=True, inplace=True)

    df0_ = pd.concat({"line0": df0_})

    # join all
    ddf = pd.concat([df0_, df_])

    ddf["z"] = ddf.z.values.astype(int)
    ddf["tag"] = ddf.tag.values.astype(int)

    # check orientation
    r0 = ddf.loc["line0"]

    if not shapely.geometry.LinearRing(r0[["lon", "lat"]].values).is_ccw:

        rf0 = ddf.loc["line0"].iloc[::-1].reset_index(drop=True)
        ddf.loc["line0"] = rf0.values

    ddf = ddf.apply(pd.to_numeric)

    return ddf


def make_bgmesh(dem, res_min, res_max):

    # scale bathymetry
    try:
        b = dem.adjusted.to_dataframe()
    except:
        b = dem.elevation.to_dataframe()

    b.columns = ["z"]

    b[b.z >= -10] = -1.0e-4  # normalize to only negative values

    b.z = np.sqrt(-b.z) / 0.5  # scale

    # adjust scale

    bg = b.z.values

    a2 = (bg - bg.min()) / (bg.max() - bg.min())

    d2 = res_min + a2 * (res_max - res_min)

    b["d2"] = d2

    nodes = b.reset_index()

    nodes["z"] = 0

    x = dem.longitude.values
    y = dem.latitude.values

    quad = MakeFacesVectorized(y.shape[0], x.shape[0])  # get element structure from array
    elems = pd.DataFrame(quad, columns=["a", "b", "c", "d"])

    df = to_df(elems, nodes)

    return df


def make_bgmesh_gradient(dem, res_min, res_max):

    # scale bathymetry
    try:
        b = dem.adjusted.to_dataframe()
    except:
        b = dem.elevation.to_dataframe()

    b.columns = ["z"]

    b[b.z >= -10] = -1.0e-4  # normalize to only negative values

    b.z = np.sqrt(-b.z) / 0.5  # scale

    # compute gradient

    sdem = b.z.values.reshape(dem.adjusted.shape)

    vgrad = np.gradient(sdem)

    mag = np.sqrt(vgrad[0] ** 2 + vgrad[1] ** 2)

    # adjust scale

    bg = mag.flatten()

    a2 = (bg - bg.min()) / (bg.max() - bg.min())

    d2 = res_min + a2 * (res_max - res_min)

    b["d2"] = d2

    nodes = b.reset_index()

    nodes["z"] = 0

    x = dem.longitude.values
    y = dem.latitude.values

    quad = MakeFacesVectorized(y.shape[0], x.shape[0])  # get element structure from array
    elems = pd.DataFrame(quad, columns=["a", "b", "c", "d"])

    df = to_df(elems, nodes)

    return df


def make_gmsh(df, **kwargs):

    logger.info("create grid")

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()
    model.add("schism")

    #    gmsh.option.setNumber("General.Terminal", 1)

    interpolate = kwargs.get("interpolate", False)
    if interpolate:
        df = gset(df, **kwargs)
    else:
        df = df
    lc = kwargs.get("lc", 0.5)

    df["lc"] = lc

    # save boundary configuration for Line0
    df_ = df.loc[df.tag != "island"].reset_index(drop=True)  # all external contours

    if not df_.empty:

        # store xy in a DataFrame
        dic = {}
        for k, d in df_.iterrows():
            out = pd.DataFrame(d.geometry.coords[:], columns=["lon", "lat"])
            out["lindex"] = d.lindex
            out = out.drop_duplicates(["lon", "lat"])
            if d.tag == "land":  # drop end points in favor or open tag
                out = out[1:-1]
            dic.update({"line{}".format(k): out})
        o1 = pd.concat(dic, axis=0).droplevel(0).reset_index(drop=True)
        o1 = o1.drop_duplicates(["lon", "lat"])

        # Do linemerge of outer contours
        lss = df_.geometry.values
        merged = shapely.ops.linemerge(lss)
        o2 = pd.DataFrame({"lon": merged.xy[0], "lat": merged.xy[1]})  # convert to DataFrame
        o2 = o2.drop_duplicates()

        rb0 = o2.merge(o1)  # merge to transfer the lindex

        rb0["z"] = 0

        # check orientation
        if not shapely.geometry.LinearRing(rb0[["lon", "lat"]].values).is_ccw:

            rb0_ = rb0.iloc[::-1].reset_index(drop=True)
            rb0 = rb0_

        rb0["lc"] = lc

        logger.info("Compute edges")
        edges = [
            list(a) for a in zip(np.arange(rb0.shape[0]), np.arange(rb0.shape[0]) + 1, rb0.lindex.values)
        ]  # outer
        edges[-1][1] = 0

        # sort (duplicated bounds)
        edges = pd.DataFrame(edges, columns=["index", "ie", "lindex"])

        edges["lindex1"] = (
            pd.concat([rb0.loc[1:, "lindex"], rb0.loc[0:0, "lindex"]]).reset_index(drop=True).astype(int)
        )
        edges["que"] = np.where(
            ((edges["lindex"] != edges["lindex1"]) & (edges["lindex"] > 0) & (edges["lindex"] < 1000)),
            edges["lindex1"],
            edges["lindex"],
        )

        edges = edges.reset_index().loc[:, ["index", "ie", "que"]]

        rb0["bounds"] = edges.loc[:, ["index", "ie"]].values.tolist()

        # get boundary types
        land_lines = {
            your_key: edges.loc[edges.que == your_key].index.values
            for your_key in [x for x in edges.que.unique() if x > 1000]
        }
        open_lines = {
            your_key: edges.loc[edges.que == your_key].index.values
            for your_key in [x for x in edges.que.unique() if x < 1000]
        }

        # join
        ols = [j for i in list(open_lines.values()) for j in i]
        lls = [j for i in list(land_lines.values()) for j in i]

    logger.info("Define geometry")

    loops = []
    islands = []
    all_lines = []

    ltag = 1

    for row in rb0.itertuples(index=True, name="Pandas"):
        factory.addPoint(
            getattr(row, "lon"), getattr(row, "lat"), getattr(row, "z"), getattr(row, "lc"), getattr(row, "Index")
        )
    for row in rb0.itertuples(index=True, name="Pandas"):
        factory.addLine(getattr(row, "bounds")[0], getattr(row, "bounds")[1], getattr(row, "Index"))

    lines = rb0.index.values
    all_lines.append(lines)

    tag = rb0.index.values[-1]

    factory.addCurveLoop(lines, tag=ltag)
    # print(loop)
    loops.append(ltag)
    all_lines.append(lines)

    tag += 1
    ltag += 1

    for k, d in df.loc[df.tag == "island"].iterrows():

        rb = pd.DataFrame(d.geometry.coords[:], columns=["lon", "lat"])
        rb["z"] = 0
        rb["lc"] = lc
        rb = rb.drop_duplicates(["lon", "lat"])

        if not shapely.geometry.LinearRing(rb[["lon", "lat"]].values).is_ccw:  # check for clockwise orientation
            rb_ = rb.iloc[::-1].reset_index(drop=True)
            rb = rb_

        rb.index = rb.index + tag
        rb["bounds"] = [[i, i + 1] for i in rb.index]
        rb["bounds"] = rb.bounds.values.tolist()[:-1] + [[rb.index[-1], rb.index[0]]]  # fix last one

        for row in rb.itertuples(index=True, name="Pandas"):
            factory.addPoint(
                getattr(row, "lon"), getattr(row, "lat"), getattr(row, "z"), getattr(row, "lc"), getattr(row, "Index")
            )
        for row in rb.itertuples(index=True, name="Pandas"):
            factory.addLine(getattr(row, "bounds")[0], getattr(row, "bounds")[1], getattr(row, "Index"))

        lines = rb.index.values
        all_lines.append(lines)

        tag = rb.index.values[-1] + 1

        factory.addCurveLoop(lines, tag=ltag)
        #    print(tag)
        loops.append(ltag)

        islands.append(lines)
        all_lines.append(lines)

        tag += 1
        ltag += 1

    factory.addPlaneSurface(loops)
    logger.info("synchronize")
    factory.synchronize()

    ## Group open boundaries lines
    for key, values in open_lines.items():
        gmsh.model.addPhysicalGroup(1, values, 1000 - int(key))

    ## Group land boundaries lines
    for key, values in land_lines.items():
        gmsh.model.addPhysicalGroup(1, values, int(key))

    ntag = 1
    for k in tqdm(range(len(islands))):
        gmsh.model.addPhysicalGroup(1, islands[k], 2000 + ntag)
        ntag += 1

    ps = gmsh.model.addPhysicalGroup(2, [1])
    gmsh.model.setPhysicalName(2, ps, "MyMesh")

    flat_list = [item for sublist in all_lines for item in sublist]
    ols = [j for i in list(open_lines.values()) for j in i]
    lists = [x for x in flat_list if x not in ols]

    model.mesh.field.add("Distance", 1)
    model.mesh.field.setNumbers(1, "CurvesList", lists)

    SizeMin = kwargs.get("SizeMin", 0.1)
    SizeMax = kwargs.get("SizeMax", 0.5)
    DistMin = kwargs.get("DistMin", 0.01)
    DistMax = kwargs.get("DistMax", 0.2)

    model.mesh.field.add("Threshold", 2)
    model.mesh.field.setNumber(2, "InField", 1)
    model.mesh.field.setNumber(2, "SizeMin", SizeMin)
    model.mesh.field.setNumber(2, "SizeMax", SizeMax)
    model.mesh.field.setNumber(2, "DistMin", DistMin)
    model.mesh.field.setNumber(2, "DistMax", DistMax)

    # Set bgmesh
    bgmesh = kwargs.get("bgmesh", None)

    if bgmesh:

        if bgmesh == "auto":

            try:

                logger.info("Read DEM")

                lon_min = df.bounds.minx.min()
                lon_max = df.bounds.maxx.max()
                lat_min = df.bounds.miny.min()
                lat_max = df.bounds.maxy.max()

                dem = pdem.dem(lon_min=lon_min, lon_max=lon_max, lat_min=lat_min, lat_max=lat_max, **kwargs)

                res_min = kwargs.get("resolution_min", 0.01)
                res_max = kwargs.get("resolution_max", 0.5)

                logger.info("Evaluate bgmesh")
                w = make_bgmesh(dem.Dataset, res_min, res_max)

                path = kwargs.get("rpath", ".")

                if not os.path.exists(path):  # check if run folder exists
                    os.makedirs(path)

                logger.info("Save bgmesh to {}/bgmesh/bgmesh.pos".format(path))

                fpos = path + "/bgmesh/bgmesh.pos"
                to_sq(w, fpos)  # save bgmesh

                kwargs.update({"bgmesh": fpos})

                model.mesh.field.setNumber(2, "StopAtDistMax", 1)

                # Merge a post-processing view containing the target anisotropic mesh sizes
                gmsh.merge(fpos)

                model.mesh.field.add("PostView", 3)
                model.mesh.field.setNumber(3, "ViewIndex", 0)

                model.mesh.field.add("Min", 4)
                model.mesh.field.setNumbers(4, "FieldsList", [2, 3])

                model.mesh.field.setAsBackgroundMesh(4)

            except:

                logger.warning("bgmesh failed... continuing without background mesh size")

                model.mesh.field.setAsBackgroundMesh(2)

        elif bgmesh.endswith(".pos"):

            gmsh.merge(bgmesh)

            model.mesh.field.setNumber(2, "StopAtDistMax", 1)

            model.mesh.field.add("PostView", 3)
            model.mesh.field.setNumber(3, "ViewIndex", 0)

            model.mesh.field.add("Min", 4)
            model.mesh.field.setNumbers(4, "FieldsList", [2, 3])

            model.mesh.field.setAsBackgroundMesh(4)

    else:

        model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    logger.info("execute gmsh")

    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    rpath = kwargs.get("rpath", ".")

    logger.info("save mesh")
    #    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(rpath + "/gmsh/mymesh.msh")

    #    gmsh.write('mymesh.vtk')

    gmsh.finalize()
