import numpy as np
import pandas as pd
import pyresample
import xarray as xr
import os
from tqdm.auto import tqdm

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
from pyposeidon.utils.fix import dem_range, resample
import logging

logger = logging.getLogger(__name__)


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
    ci = 4 * R**2 / (ui**2 + vi**2 + 4 * R**2)

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

    m = mesh.Dataset
    m = m.assign({"x": m["SCHISM_hgrid_node_x"], "y": m["SCHISM_hgrid_node_y"]})
    # Stereo -> lat/lon
    clon, clat = to_lat_lon(x0, y0)
    m["x"].data = clon
    m["y"].data = clat

    # Select DEM
    #    try:
    #        dm = dem.adjusted.to_dataframe()
    #    except:
    #        dm = dem.elevation.to_dataframe()

    #    lon = dem.longitude.values
    #    lat = dem.latitude.values

    #    X, Y = np.meshgrid(lon, lat)
    # Stereo -> lat/lon
    #    clon, clat = to_lat_lon(x0, y0)

    # resample bathymetry
    #    gdem = dm.values.flatten()

    #    orig = pyresample.geometry.SwathDefinition(lons=X.flatten(), lats=Y.flatten())  # original bathymetry points
    #    targ = pyresample.geometry.SwathDefinition(lons=clon, lats=clat)  # wet points

    #    bw = pyresample.kd_tree.resample_nearest(orig, gdem, targ, radius_of_influence=50000, fill_value=0)

    dem_on_mesh(m, dem)

    bw = m.depth.data

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


def fillv(dem, perms, m, buffer=0.0):

    for (i1, i2), (j1, j2) in tqdm(perms, total=len(perms)):

        lon1 = dem.longitude.data[i1:i2][0]
        lon2 = dem.longitude.data[i1:i2][-1]
        lat1 = dem.latitude.data[j1:j2][0]
        lat2 = dem.latitude.data[j1:j2][-1]

        # buffer lat/lon
        blon1 = lon1 - buffer
        blon2 = lon2 + buffer
        blat1 = lat1 - buffer
        blat2 = lat2 + buffer

        #    de = dem.sel(lon=slice(blon1,blon2)).sel(lat=slice(blat1,blat2))
        de = dem_range(dem, blon1, blon2, blat1, blat2)

        # subset mesh
        indices_of_nodes_in_bbox = np.where(
            (m.y >= lat1 - buffer / 2)
            & (m.y <= lat2 + buffer / 2)
            & (m.x >= lon1 - buffer / 2)
            & (m.x <= lon2 + buffer / 2)
        )[0]

        bm = m.isel(nSCHISM_hgrid_node=indices_of_nodes_in_bbox)
        ids = np.argwhere(np.isnan(bm.depth.values)).flatten()

        grid_x, grid_y = bm.x.data, bm.y.data

        bd = resample(de, grid_x, grid_y, var="adjusted", wet=True, flag=0, function="gauss")

        m["depth"].loc[dict(nSCHISM_hgrid_node=indices_of_nodes_in_bbox)] = -bd


def dem_on_mesh(mesh, dem):

    ilats = dem.elevation.chunk("auto").chunks[0]
    ilons = dem.elevation.chunk("auto").chunks[1]

    if len(ilons) == 1:
        ilons = (int(ilons[0] / 2), int(ilons[0] / 2))

    idx = [sum(ilons[:i]) for i in range(len(ilons) + 1)]
    jdx = [sum(ilats[:i]) for i in range(len(ilats) + 1)]

    blon = list(zip(idx[:-1], idx[1:]))
    blat = list(zip(jdx[:-1], jdx[1:]))

    perms = [(x, y) for x in blon for y in blat]

    fillv(dem, perms, mesh, buffer=5)
