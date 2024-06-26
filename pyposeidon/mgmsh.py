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
from tqdm.auto import tqdm
import sys
import subprocess
import shapely
import shlex

from pyposeidon.utils.spline import use_spline
from pyposeidon.utils.global_bgmesh import make_bgmesh_global
from pyposeidon.utils.stereo import to_lat_lon
from pyposeidon.utils.pos import to_global_pos
import pyposeidon.dem as pdem
from pyposeidon.tools import orient
from pyposeidon.dem_tools import compute_dem_gradient, make_bgmesh

import multiprocessing

NCORES = max(1, min(multiprocessing.cpu_count() - 1, 20))

from joblib import Parallel, delayed, parallel_backend

import logging

logger = logging.getLogger(__name__)


def get_ibounds(df, mm):
    if df.shape[0] == 0:
        disable = True
    else:
        disable = False

    ibounds = []

    for row in tqdm(df.itertuples(index=True, name="Pandas"), total=df.shape[0], disable=disable):
        inodes, xyz = mm.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))
        db = pd.DataFrame({"node": inodes - 1})
        db["type"] = "island"
        db["id"] = -(getattr(row, "Index") + 1)

        ibounds.append(db)

    return ibounds


def read_msh(filename, **kwargs):
    import gmsh

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()

    gmsh.open(filename)

    logger.info("Analyzing mesh")

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
    obs = bgs.loc[bgs.tag < 1000]

    logger.info("No of Open boundaries: {}".format(obs.shape[0]))

    for row in obs.itertuples(index=True, name="Pandas"):
        onodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))

        db = pd.DataFrame({"node": onodes - 1})
        db["type"] = "open"
        db["id"] = getattr(row, "Index") + 1

        bounds.append(db)

    # land boundaries type

    lbs = bgs.loc[(bgs.tag > 1000) & (bgs.tag < 2000)]
    lbs.reset_index(inplace=True, drop=True)

    logger.info("No of Land boundaries: {}".format(lbs.shape[0]))

    for row in lbs.itertuples(index=True, name="Pandas"):
        lnodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"), tag=getattr(row, "tag"))

        db = pd.DataFrame({"node": lnodes - 1})
        db["type"] = "land"
        db["id"] = 1000 + (getattr(row, "Index") + 1)

        bounds.append(db)

    # islands

    ibs = bgs.loc[bgs.tag > 2000]
    ibs.reset_index(inplace=True, drop=True)

    logger.info("No. of Islands: {}".format(ibs.shape[0]))

    # ===============================
    # parallize
    # ===============================
    mm = gmsh.model.mesh
    num_partitions = NCORES if ibs.shape[0] > NCORES else 1  # number of partitions to split dataframe
    df_split = np.array_split(ibs, num_partitions)

    with parallel_backend("threading", n_jobs=num_partitions):
        results = Parallel()(delayed(get_ibounds)(_, mm) for _ in df_split)

    bounds.extend([j for i in results for j in i])

    # ===============================

    if bounds != []:
        bnodes = pd.concat(bounds).reset_index(drop=True)

        bnodes.index.name = "bnodes"

        bnodes = bnodes.drop_duplicates("node")

        bnodes["id"] = bnodes.id.astype(int)

    else:
        bnodes = pd.DataFrame({})

    gglobal = kwargs.get("gglobal", False)
    if gglobal:
        use_bindings = kwargs.get("use_bindings", True)
        if not use_bindings:
            bnodes["type"] = "island"  # Fix for binary run and GLOBAL. CHECK
            bnodes["id"] = [-i if i > 0 else i for i in bnodes.id.values]

        bnodes = bnodes.sort_values(["type", "id", "node"]).reset_index(drop=True)  # sort
        bnodes.index.name = "bnodes"

        # check orientation
        # if use_bindings:
        # nodes, tria = orient(nodes, tria, x="x", y="y")
        #        else:
        bgmesh = kwargs.get("bgmesh", None)
        if not bgmesh:
            tria = tria.reindex(columns=["a", "c", "b"])

    # check if global and reproject
    if gglobal:  # convert to lat/lon
        if nodes.z.any() != 0:
            xd, yd = to_lat_lon(nodes.x, nodes.y, nodes.z)
            nodes["x"] = xd
            nodes["y"] = yd
        else:
            xd, yd = to_lat_lon(nodes.x, nodes.y)
            nodes["x"] = xd
            nodes["y"] = yd

    mesh = pd.DataFrame({"lon": nodes.x, "lat": nodes.y})

    tri3 = tria.values

    logger.info("Finalizing Dataset")

    ## make dataset
    els = xr.DataArray(
        tri3,
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        mesh.loc[:, ["lon", "lat"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "lon": "SCHISM_hgrid_node_x",
                "lat": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    dep = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    gr = xr.merge([nod, els, dep, bnodes.to_xarray()])

    gmsh.finalize()

    return gr


def get(contours, **kwargs):
    """
    Create a `gmsh` mesh.

    !!! danger ""
        Due to a limitation of the Library rendering the docstrings, all arguments are marked
        as `required`, nevertheless they are all `Optional`.

    Args:
        contours GeoDataFrame: Provide boundaries and metadata.
        rpath str: Path for output. Defaults to `"."`.
        use_bindings bool: Flag for using python API as opposed to binary. Defaults to `True`.
        dem_source str: Path or url to bathymetric data.
        bgmesh str: Path to a mesh scale file. Defaults to `None`.
        setup_only bool: Flag for setup only (no execution). Defaults to `False`.
    """

    logger.info("Creating mesh with GMSH\n")

    rpath = kwargs.get("rpath", ".")

    if not os.path.exists(rpath):
        os.makedirs(rpath)

    gpath = os.path.join(rpath, "gmsh")
    if not os.path.exists(gpath):
        os.makedirs(gpath)

    use_bindings = kwargs.get("use_bindings", True)
    setup_only = kwargs.get("setup_only", False)
    bgmesh = kwargs.get("bgmesh", None)

    if bgmesh is None:
        dem_source = kwargs.get("dem_source", None)
        if dem_source:
            bgmesh = "auto"
            kwargs.update({"bgmesh": "auto"})

    gglobal = kwargs.get("gglobal", False)

    dh1 = None
    dh2 = None

    if bgmesh == "auto":

        bg_dem = kwargs.get("bg_dem", True)

        if bg_dem:
            try:
                rpath = kwargs.get("rpath", ".")

                if not os.path.exists(rpath + "/gmsh/"):  # check if run folder exists
                    os.makedirs(rpath + "/gmsh/")

                fpos = rpath + "/gmsh/bgmesh1.pos"

                if gglobal:
                    dem = pdem.Dem(**kwargs)

                    nds, lms = make_bgmesh_global(contours, fpos, dem.Dataset, **kwargs)
                    dh1 = to_global_pos(nds, lms, fpos, **kwargs)

                else:
                    dh1 = make_bgmesh(contours, fpos, **kwargs)

                kwargs.update({"bgmesh1": fpos})

            except OSError as e:
                logger.warning("bgmesh1 failed... continuing without background mesh size")
                dh1 = None
                kwargs.update({"bgmesh1": None})
        else:
            dh1 = None

        # same for dem gradient
        bg_dem_grad = kwargs.get("bg_dem_grad", True)

        if bg_dem_grad:
            try:
                fpos = rpath + "/gmsh/bgmesh2.pos"

                slope_parameter = kwargs.get("slope_parameter", 20)
                filter_quotient = kwargs.get("filter_quotient", 50)
                min_edge_length = kwargs.get("min_edge_length", 0.3)
                max_edge_length = kwargs.get("max_edge_length", 2)
                min_elevation_cutoff = kwargs.get("min_elevation_cutoff", -50.0)

                dem = pdem.Dem(**kwargs)
                gdem = compute_dem_gradient(
                    dem.Dataset,
                    slope_parameter=slope_parameter,
                    filter_quotient=filter_quotient,
                    min_edge_length=min_edge_length,
                    max_edge_length=max_edge_length,
                    min_elevation_cutoff=-50.0,
                )

                kwargs.update({"var": "gradient"})

                if gglobal:
                    nds, lms = make_bgmesh_global(contours, fpos, gdem, scale=False, **kwargs)
                    dh2 = to_global_pos(nds, lms, fpos, **kwargs)
                else:
                    dh2 = make_bgmesh(contours, fpos, gdem, scale=False, **kwargs)

                kwargs.update({"bgmesh2": fpos})

            except OSError as e:
                logger.warning("bgmesh2 failed... continuing without background mesh size")
                dh2 = None
                kwargs.update({"bgmesh2": None})

        else:
            dh2 = None

    if use_bindings:
        logger.info("Using python bindings")
        if gglobal:
            gr = make_gmsh_3d(contours, **kwargs)
        else:
            make_gmsh(contours, **kwargs)
            gr = read_msh(rpath + "/gmsh/mymesh.msh", **kwargs)

    else:
        to_geo(contours, **kwargs)
        if not setup_only:
            logger.info("Using GMSH binary")
            gmsh_execute(**kwargs)
            gr = read_msh(rpath + "/gmsh/mymesh.msh", **kwargs)
        else:
            gr = None

    bg = [dh1, dh2]
    return gr, bg


def gset(df, **kwargs):
    logger.info("Interpolating coastal points")

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
        di = use_spline(curve, ds=0.01, method="slinear")
        di["z"] = df.loc[contour].z.values[0]
        di["tag"] = df.loc[contour].tag.values[0].astype(int)
        di["lc"] = df.loc[contour].lc.values[0]
        ndfsfs.update({contour: di.drop_duplicates(["lon", "lat"])})

    df_ = pd.concat(ndfsfs, axis=0)
    df_["z"] = df_.z.values.astype(int)
    df_["tag"] = df_.tag.values.astype(int)

    # Line0

    logger.info("Setting outermost boundary")

    df0 = df.loc["line0"]

    mtag = df0.tag.min()
    mtag = mtag.astype(int)

    nd0 = {}
    for ic in tqdm(range(mtag, 0)):
        contour = df0.tag == ic
        curve = df0.loc[contour, ["lon", "lat"]].reset_index(drop=True)
        #    curve = pd.concat([curve,curve.loc[0:0]]).reset_index(drop=True)
        di = use_spline(curve, ds=0.01, method="slinear")
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


def to_geo(df, **kwargs):
    gglobal = kwargs.get("gglobal", False)

    lc = kwargs.get("lc", 0.5)
    df["lc"] = lc

    bspline = kwargs.get("bspline", False)

    # save boundary configuration for Line0
    df_ = df.loc[df.tag != "island"].reset_index(drop=True)  # all external contours

    if not df_.empty:
        rb0, land_lines, open_lines = outer_boundary(df, **kwargs)
    else:
        rb0 = pd.DataFrame({})
        open_lines = {}
        land_lines = {}

    ptag = 0
    ptag0 = 1
    ltag = 0
    loops = []

    # ... and save it to disk
    rpath = kwargs.get("rpath", ".")

    logger.info("Saving geo file")
    filename = rpath + "/gmsh/mymesh.geo"

    with open(filename, "w") as f:
        if gglobal:
            f.write("Point(1) = {{0,0,0,{}}};\n".format(lc))
            f.write("Point(2) = {{1,0,0,{}}};\n".format(lc))
            f.write("PolarSphere (1) = {1,2};\n")
            ptag = 2
            ptag0 = 3
            ltag = 1
            loops = []

        # outer contour
        points = []
        for row in rb0.itertuples(index=True, name="Pandas"):
            ptag += 1
            f.write(
                "Point({}) = {{{},{},{},{}}};\n".format(
                    ptag,
                    getattr(row, "lon"),
                    getattr(row, "lat"),
                    getattr(row, "z"),
                    getattr(row, "lc"),
                )
            )
            points.append(ptag)

        if points:
            points = points + [points[0]]

            if not bspline:
                t2 = [[points[l], points[l + 1]] for l in range(0, len(points) - 1)]
                lines = []
                for [a, b] in t2:
                    ltag += 1
                    f.write("Line({}) = {{{}, {}}};\n".format(ltag, a, b))
                    lines.append(ltag)
                ltag += 1

            else:
                ltag += 1
                f.write("BSpline ({}) = {}".format(ltag, "{"))
                f.write(",".join(map(str, points)))
                f.write("{};\n".format("}"))
                lines = [ltag]
                ltag += 1

            f.write("Line Loop ({}) = {}".format(ltag, "{"))
            f.write(",".join(map(str, lines)))
            f.write("{};\n".format("}"))

            loops.append(ltag)

        if not bspline:
            ## Group open boundaries lines
            for key, values in open_lines.items():
                f.write("Physical Line  ({}) = {}".format(key, "{"))
                f.write(",".join(map(str, values + 1)))
                f.write("{};\n".format("}"))

            ## Group land boundaries lines
            for key, values in land_lines.items():
                f.write("Physical Line  ({}) = {}".format(key, "{"))
                f.write(",".join(map(str, values + 1)))
                f.write("{};\n".format("}"))

        # The rest (islands)
        pl = 2000

        dfi = df.loc[df.tag == "island"]

        for k, d in dfi.iterrows():
            points = []

            rb = pd.DataFrame(d.geometry.coords[:], columns=["lon", "lat"])
            rb["z"] = 0
            rb["lc"] = lc
            rb = rb.drop_duplicates(["lon", "lat"])

            if not shapely.geometry.LinearRing(rb[["lon", "lat"]].values).is_ccw:  # check for clockwise orientation
                rb_ = rb.iloc[::-1].reset_index(drop=True)
                rb = rb_

            for row in rb.itertuples(index=True, name="Pandas"):
                ptag += 1
                f.write(
                    "Point({}) = {{{},{},{},{}}};\n".format(
                        ptag,
                        getattr(row, "lon"),
                        getattr(row, "lat"),
                        getattr(row, "z"),
                        getattr(row, "lc"),
                    )
                )
                points.append(ptag)

            points = points + [points[0]]

            if not bspline:
                t2 = [[points[l], points[l + 1]] for l in range(0, len(points) - 1)]
                lines = []
                for [a, b] in t2:
                    ltag += 1
                    f.write("Line({}) = {{{}, {}}};\n".format(ltag, a, b))
                    lines.append(ltag)
                ltag += 1

            else:
                ltag += 1
                f.write("BSpline ({}) = {}".format(ltag, "{"))
                f.write(",".join(map(str, points)))
                f.write("{};\n".format("}"))
                lines = [ltag]
                ltag += 1

            f.write("Line Loop ({}) = {}".format(ltag, "{"))
            f.write(",".join(map(str, lines)))
            f.write("{};\n".format("}"))

            loops.append(ltag)

            if not bspline:
                pl += 1
                #            f.write("Physical Line ({}) = {{{}}};\n".format(pl, ltag))
                f.write("Physical Line  ({}) = {}".format(pl, "{"))
                f.write(",".join(map(str, lines)))
                f.write("{};\n".format("}"))

        f.write("Plane Surface (2) = {")
        f.write(",".join(map(str, loops)))
        f.write("{};\n".format("}"))

        f.write("Physical Surface(1) = {2};\n")
        f.write("Field[1] = Distance;\n")

        if not gglobal:
            curve_list = [x for x in np.arange(1, ltag + 1) if x not in loops]
        else:
            curve_list = [x for x in np.arange(2, ltag + 1) if x not in loops]

        if not bspline:
            f.write("Field[1].CurvesList = {{{}}};\n".format(",".join(map(str, curve_list))))
        else:
            f.write("Field[1].NodesList =\n{")
            f.write(",".join(map(str, np.arange(ptag0, ptag + 1))))
            f.write("};\n")

        SizeMin = kwargs.get("SizeMin", 0.1)
        SizeMax = kwargs.get("SizeMax", 0.5)
        DistMin = kwargs.get("DistMin", 0.01)
        DistMax = kwargs.get("DistMax", 0.2)

        f.write("Field[2] = Threshold;\n")
        f.write("Field[2].IField = 1;\n")
        f.write("Field[2].DistMin = {};\n".format(DistMin))
        f.write("Field[2].DistMax = {};\n".format(DistMax))
        f.write("Field[2].SizeMin = {};\n".format(SizeMin))
        f.write("Field[2].SizeMax = {};\n".format(SizeMax))

        f.write("Field[2].StopAtDistMax = 1;\n")

        bgmesh = kwargs.get("bgmesh", None)
        if bgmesh == "auto":
            bgmesh = None
        bgmesh1 = kwargs.get("bgmesh1", None)
        bgmesh2 = kwargs.get("bgmesh2", None)

        if any([bgmesh, bgmesh1, bgmesh2]):

            ib = 2
            vi = -1

            for bgm in [bgmesh, bgmesh1, bgmesh2]:

                if bgm:

                    ib += 1
                    vi += 1
                    # get full path
                    bgmesh_path = os.path.abspath(bgm)

                    f.write('Merge "{}";\n'.format(bgmesh_path))

                    f.write(f"Field[{ib}] = PostView;\n")
                    f.write(f"Field[{ib}].ViewIndex = {vi};\n")

            flist = np.arange(2, ib + 1)
            f.write(f"Field[{ib+1}] = Min;\n")
            f.write(f"Field[{ib+1}].FieldsList = {{{', '.join(map(str,flist))}}};\n")

            f.write(f"Background Field = {ib+1};\n")

        else:
            f.write("Background Field = 2;\n")

        MeshSizeMin = kwargs.get("MeshSizeMin", SizeMin)
        MeshSizeMax = kwargs.get("MeshSizeMax", SizeMax)

        f.write(f"Mesh.MeshSizeMin= {MeshSizeMin};\n")
        f.write(f"Mesh.MeshSizeMax= {MeshSizeMax};\n")

        f.write("Mesh.MeshSizeFromPoints = 0;\n")
        f.write("Mesh.MeshSizeFromCurvature = 0;\n")
        f.write("Mesh.MeshSizeExtendFromBoundary = 0;\n")

    return


def gmsh_execute(**kwargs):
    rpath = kwargs.get("rpath", ".")
    setup_only = kwargs.get("setup_only", False)
    calc_dir = rpath + "/gmsh/"
    ncores = kwargs.get("ncores", NCORES)

    try:
        bin_path = os.environ["GMSH"]
    except:
        bin_path = kwargs.get("gpath", None)

    if bin_path is None:
        # ------------------------------------------------------------------------------
        logger.warning("gmsh executable path (gpath) not given -> using default \n")
        # ------------------------------------------------------------------------------
        bin_path = "gmsh"

    if not setup_only:
        # ---------------------------------
        logger.info("Executing gmsh\n")
        # ---------------------------------

        gglobal = kwargs.get("gglobal", False)
        if gglobal:
            dim = -3
        else:
            dim = -2

        gmsh_args = kwargs.get("gmsh_args", {})

        gargs = ""
        for k, v in list(gmsh_args.items()):
            if k[0] != "-":
                k = "-" + k
            gargs = gargs + " ".join([k, str(v)]) + " "

        # execute gmsh
        cmd = "{} {} {} {}".format(bin_path, gargs, dim, "mymesh.geo")

        proc = subprocess.run(
            shlex.split(cmd),
            check=False,
            capture_output=True,
            text=True,
            cwd=calc_dir,
            # bufsize=1,
        )

        with open(os.path.join(calc_dir, "myerr.log"), "w") as fd:
            fd.write(proc.stderr)
        with open(os.path.join(calc_dir, "myrun.log"), "w") as fd:
            fd.write(proc.stdout)

        if proc.returncode or "Error" in proc.stderr:
            # ---------------------------------------------------------------------
            logger.error("gmsh failed to execute\n")
            # ---------------------------------------------------------------------
            proc.check_returncode()
            raise subprocess.CalledProcessError(
                cmd=cmd, output=proc.stdout, stderr=proc.stderr, returncode=proc.returncode
            )
        else:
            logger.info("gmsh finished successfully\n")

    return


def outer_boundary(df, **kwargs):
    lc = kwargs.get("lc", 0.5)

    df_ = df.loc[df.tag != "island"].reset_index(drop=True)  # all external contours

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
    merged = shapely.ops.linemerge(list(lss))
    o2 = pd.DataFrame({"lon": merged.xy[0], "lat": merged.xy[1]})  # convert to DataFrame
    o2 = o2.drop_duplicates()

    rb0 = o2.merge(o1)  # merge to transfer the lindex

    rb0["z"] = 0

    # check orientation
    if not shapely.geometry.LinearRing(rb0[["lon", "lat"]].values).is_ccw:
        rb0_ = rb0.iloc[::-1].reset_index(drop=True)
        rb0 = rb0_

    rb0["lc"] = lc

    logger.info("Computing edges")
    edges = [list(a) for a in zip(np.arange(rb0.shape[0]), np.arange(rb0.shape[0]) + 1, rb0.lindex.values)]  # outer
    edges[-1][1] = 0

    # sort (duplicated bounds)
    edges = pd.DataFrame(edges, columns=["index", "ie", "lindex"])

    edges["lindex1"] = pd.concat([rb0.loc[1:, "lindex"], rb0.loc[0:0, "lindex"]]).reset_index(drop=True).astype(int)
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

    return rb0, land_lines, open_lines


def make_gmsh(df, **kwargs):
    import gmsh

    logger.info("Creating mesh")

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
        rb0, land_lines, open_lines = outer_boundary(df, **kwargs)
        # join
        ols = [j for i in list(open_lines.values()) for j in i]
        lls = [j for i in list(land_lines.values()) for j in i]

    else:
        rb0 = pd.DataFrame({})
        open_lines = {}
        land_lines = {}

    logger.info("Defining geometry")

    loops = []
    islands = []
    all_lines = []

    ltag = 1
    tag = 0

    if not df_.empty:
        for row in rb0.itertuples(index=True, name="Pandas"):
            factory.addPoint(
                getattr(row, "lon"),
                getattr(row, "lat"),
                getattr(row, "z"),
                getattr(row, "lc"),
                getattr(row, "Index"),
            )
        for row in rb0.itertuples(index=True, name="Pandas"):
            factory.addLine(
                getattr(row, "bounds")[0],
                getattr(row, "bounds")[1],
                getattr(row, "Index"),
            )

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
                getattr(row, "lon"),
                getattr(row, "lat"),
                getattr(row, "z"),
                getattr(row, "lc"),
                getattr(row, "Index"),
            )
        for row in rb.itertuples(index=True, name="Pandas"):
            factory.addLine(
                getattr(row, "bounds")[0],
                getattr(row, "bounds")[1],
                getattr(row, "Index"),
            )

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
    logger.info("Synchronizing")
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
    model.mesh.field.setNumber(2, "StopAtDistMax", 1)

    # Set bgmesh
    bgmesh = kwargs.get("bgmesh", None)
    if bgmesh == "auto":
        bgmesh = None
    bgmesh1 = kwargs.get("bgmesh1", None)
    bgmesh2 = kwargs.get("bgmesh2", None)

    if any([bgmesh, bgmesh1, bgmesh2]):

        ib = 2
        vi = -1

        for bgm in [bgmesh, bgmesh1, bgmesh2]:

            if bgm:

                ib += 1
                vi += 1

                gmsh.merge(bgm)

                model.mesh.field.add("PostView", ib)
                model.mesh.field.setNumber(ib, "ViewIndex", vi)

        blist = np.arange(2, ib + 1).tolist()
        model.mesh.field.add("Min", ib + 1)
        model.mesh.field.setNumbers(ib + 1, "FieldsList", blist)

        model.mesh.field.setAsBackgroundMesh(ib + 1)

    else:
        model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    MeshSizeMin = kwargs.get("MeshSizeMin", SizeMin)
    MeshSizeMax = kwargs.get("MeshSizeMax", SizeMax)

    if MeshSizeMin is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMin", MeshSizeMin)
    if MeshSizeMax is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMax", MeshSizeMax)

    logger.info("Executing gmsh")

    gmsh.option.set_number("General.Verbosity", 0)

    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    rpath = kwargs.get("rpath", ".")

    logger.info("Saving mesh")
    #    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(rpath + "/gmsh/mymesh.msh")

    gmsh.write(rpath + "/gmsh/mymesh.vtk")

    gmsh.finalize()

    return


def make_gmsh_3d(df, **kwargs):
    import gmsh

    logger.info("Creating global mesh")

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()
    model.add("schism")

    #    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Geometry.Tolerance", 1e-20)

    interpolate = kwargs.get("interpolate", False)
    if interpolate:
        df = gset(df, **kwargs)
    else:
        df = df

    logger.info("Defining geometry")

    loops = []
    all_lines = []

    curve_tag = 1
    tag = 1

    factory.addGeometry("PolarSphere", [0, 0, 0, 1], tag=tag)

    tag += 1
    curve_tag += 1

    for k, d in df.iterrows():
        rb = pd.DataFrame(d.geometry.coords[:], columns=["lon", "lat"])
        rb["z"] = 1
        rb = rb.drop_duplicates(["lon", "lat"])

        if not shapely.geometry.LinearRing(rb[["lon", "lat"]].values).is_ccw:  # check for clockwise orientation
            rb_ = rb.iloc[::-1].reset_index(drop=True)
            rb = rb_

        rb.index = rb.index + tag
        rb["bounds"] = [[i, i + 1] for i in rb.index]
        rb["bounds"] = rb.bounds.values.tolist()[:-1] + [[rb.index[-1], rb.index[0]]]  # fix last one

        for row in rb.itertuples(index=True, name="Pandas"):
            factory.addPointOnGeometry(
                getattr(row, "z"),
                getattr(row, "lon"),
                getattr(row, "lat"),
                tag=getattr(row, "Index"),
            )

        for row in rb.itertuples(index=True, name="Pandas"):
            factory.addLine(
                getattr(row, "bounds")[0],
                getattr(row, "bounds")[1],
                getattr(row, "Index"),
            )

        lines = rb.index.values
        all_lines.append(lines)

        tag = rb.index.values[-1] + 1

        factory.addCurveLoop(lines, tag=curve_tag)
        #    print(tag)
        loops.append(curve_tag)

        tag += 1
        curve_tag += 1

    factory.addPlaneSurface(loops)
    logger.info("Synchronizing")
    factory.synchronize()

    ## Group boundaries lines
    ntag = 1
    for k in tqdm(range(len(all_lines))):
        gmsh.model.addPhysicalGroup(1, all_lines[k], 2000 + ntag)
        ntag += 1

    ps = gmsh.model.addPhysicalGroup(2, [1])
    gmsh.model.setPhysicalName(2, ps, "MyMesh")

    model.mesh.field.add("Distance", 1)
    model.mesh.field.setNumbers(1, "CurvesList", [d[1] for d in gmsh.model.getEntities(1)])

    SizeMin = kwargs.get("SizeMin", 0.01)
    SizeMax = kwargs.get("SizeMax", 0.1)
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
    if bgmesh == "auto":
        bgmesh = None
    bgmesh1 = kwargs.get("bgmesh1", None)
    bgmesh2 = kwargs.get("bgmesh2", None)

    if any([bgmesh, bgmesh1, bgmesh2]):

        ib = 2
        vi = -1

        for bgm in [bgmesh, bgmesh1, bgmesh2]:

            if bgm:

                ib += 1
                vi += 1

                gmsh.merge(bgm)

                model.mesh.field.add("PostView", ib)
                model.mesh.field.setNumber(ib, "ViewIndex", vi)

        blist = np.arange(2, ib + 1).tolist()
        model.mesh.field.add("Min", ib + 1)
        model.mesh.field.setNumbers(ib + 1, "FieldsList", blist)

        model.mesh.field.setAsBackgroundMesh(ib + 1)

    else:
        model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    MeshSizeMin = kwargs.get("MeshSizeMin", None)
    MeshSizeMax = kwargs.get("MeshSizeMax", None)

    if MeshSizeMin is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMin", MeshSizeMin)
    if MeshSizeMax is not None:
        gmsh.option.setNumber("Mesh.MeshSizeMax", MeshSizeMax)

    logger.info("Executing gmsh")

    gmsh.option.set_number("General.Verbosity", 0)

    gmsh.model.mesh.generate(3)

    logger.info("Analyzing mesh")

    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    nodes = pd.DataFrame(coord.reshape(-1, 3), columns=["x", "y", "z"])

    elementTags2, nodeTags2 = gmsh.model.mesh.getElementsByType(2)

    elems = nodeTags2.reshape(-1, 3)

    tria = pd.DataFrame(elems - 1, columns=["a", "b", "c"])

    # boundaries

    bgs = gmsh.model.getPhysicalGroups()

    bgs = pd.DataFrame(bgs[:-1], columns=["dim", "tag"])

    mm = gmsh.model.mesh
    num_partitions = NCORES if bgs.shape[0] > NCORES else 1  # number of partitions to split dataframe
    df_split = np.array_split(bgs, num_partitions)

    with parallel_backend("threading", n_jobs=num_partitions):
        results = Parallel()(delayed(get_ibounds)(_, mm) for _ in df_split)

    bounds = [j for i in results for j in i]

    if bounds != []:
        bnodes = pd.concat(bounds).reset_index(drop=True)

        bnodes.index.name = "bnodes"

        bnodes = bnodes.drop_duplicates("node")

        bnodes["id"] = bnodes.id.astype(int)

    else:
        bnodes = pd.DataFrame({})

    bnodes = bnodes.sort_values(["id", "node"]).reset_index(drop=True)  # sort
    bnodes.index.name = "bnodes"

    # check orientation
    #    nodes, tria = orient(nodes, tria, x="x", y="y")

    # convert to lat/lon
    if nodes.z.any() != 0:
        xd, yd = to_lat_lon(nodes.x, nodes.y, nodes.z)
        nodes["x"] = xd
        nodes["y"] = yd
    else:
        xd, yd = to_lat_lon(nodes.x, nodes.y)
        nodes["x"] = xd
        nodes["y"] = yd

    mesh = pd.DataFrame({"lon": nodes.x, "lat": nodes.y})

    tri3 = tria.values

    logger.info("Finalizing Dataset")

    ## make dataset
    els = xr.DataArray(
        tri3,
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        mesh.loc[:, ["lon", "lat"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "lon": "SCHISM_hgrid_node_x",
                "lat": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    dep = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    gr = xr.merge([nod, els, dep, bnodes.to_xarray()])

    # ... and save it to disk
    rpath = kwargs.get("rpath", ".")

    logger.info("Saving mesh")
    #    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(rpath + "/gmsh/mymesh.msh")

    gmsh.write(rpath + "/gmsh/mymesh.vtk")

    gmsh.finalize()

    return gr
