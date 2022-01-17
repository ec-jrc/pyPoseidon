import numpy as np
import pandas as pd
import pyresample
import xarray as xr
import os

import pyposeidon.boundary as pb
from pyposeidon.utils.stereo import to_3d, to_lat_lon, stereo_to_3d
from pyposeidon.utils.pos import to_sq, to_st
from pyposeidon.utils.hfun import to_hfun_mesh, to_hfun_grid
from pyposeidon.utils.scale import scale_dem
from pyposeidon.utils.topology import (
    MakeQuadFaces,
    quads_to_df,
    tria_to_df_3d,
    tria_to_df,
)
import pyposeidon.dem as pdem
import pyposeidon.mesh as pmesh
import logging

logger = logging.getLogger("pyposeidon")


def make_bgmesh_global(dfb, fpos, dem, **kwargs):

    mesh_generator = kwargs.get("mesh_generator", None)

    bk = dfb.copy()

    dfb = dfb.loc[dfb.length == dfb.length.max()]

    out = dfb.buffer(0.2).exterior

    R = kwargs.get("R", 1.0)

    bgm_res = kwargs.get("bgm_res", 0.01)

    # create simple mesh
    npu = kwargs.get("npu", 1000)

    d1 = np.linspace(out.bounds.minx, out.bounds.maxx, npu)
    d2 = np.linspace(out.bounds.miny, out.bounds.maxy, npu)

    ui, vi = np.meshgrid(d1, d2)

    # stereo->2D scale
    ci = 4 * R ** 2 / (ui ** 2 + vi ** 2 + 4 * R ** 2)

    # set weight field
    scale = bgm_res / ci.flatten()

    # save as bgmesh
    nodes = pd.DataFrame({"longitude": ui.flatten(), "latitude": vi.flatten(), "z": 0, "d2": scale})

    # tesselation
    quad = MakeQuadFaces(npu, npu)
    elems = pd.DataFrame(quad, columns=["a", "b", "c", "d"])

    if mesh_generator == "gmsh":

        dn = quads_to_df(elems, nodes)

        logger.info("Save pre-process bgmesh to {}".format(fpos))

        to_sq(dn, fpos)  # save bgmesh

    elif mesh_generator == "jigsaw":

        dh = xr.Dataset(
            {"h": (["longitude", "latitude"], nodes.d2.values.reshape(ui.shape))},
            coords={
                "longitude": ("longitude", d1.flatten()),
                "latitude": ("latitude", d2.flatten()),
            },
        )

        logger.info("Save pre-process bgmesh to {}".format(fpos))

        to_hfun_grid(dh, fpos)  # write bgmesh file

    rpath = kwargs.get("rpath", ".")

    logger.info("Create interim global scale mesh")

    b = pb.Boundary(geometry=dfb)
    mesh = pmesh.set(
        type="tri2d",
        boundary=b,
        mesh_generator=mesh_generator,
        rpath=rpath,
        bgmesh=fpos,
    )

    x0 = mesh.Dataset.SCHISM_hgrid_node_x.values
    y0 = mesh.Dataset.SCHISM_hgrid_node_y.values
    trii0 = mesh.Dataset.SCHISM_hgrid_face_nodes.values[:, :3]

    # Select DEM
    try:
        dm = dem.adjusted.to_dataframe()
    except:
        dm = dem.elevation.to_dataframe()

    lon = dem.longitude.values
    lat = dem.latitude.values

    X, Y = np.meshgrid(lon, lat)
    # Stereo -> lat/lon
    clon, clat = to_lat_lon(x0, y0)

    # resample bathymetry
    gdem = dm.values.flatten()

    orig = pyresample.geometry.SwathDefinition(lons=X.flatten(), lats=Y.flatten())  # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=clon, lats=clat)  # wet points

    bw = pyresample.kd_tree.resample_nearest(orig, gdem, targ, radius_of_influence=50000, fill_value=0)

    bz = pd.DataFrame({"z": bw.flatten()})

    # scale
    res_min = kwargs.get("resolution_min", 0.01)
    res_max = kwargs.get("resolution_max", 0.5)

    nodes = scale_dem(bz, res_min, res_max, **kwargs)

    nodes["u"] = x0
    nodes["v"] = y0

    elems = pd.DataFrame(trii0, columns=["a", "b", "c"])

    dfb = bk

    return nodes, elems
