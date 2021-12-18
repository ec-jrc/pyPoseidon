import pandas as pd
import geopandas as gp
import shapely
import numpy as np
import xarray as xr
from .limgrad import *
import matplotlib
from pyposeidon.utils.stereo import to_lat_lon, to_stereo
from pyposeidon.utils.topology import MakeTriangleFaces, MakeTriangleFaces_periodic
from pyposeidon.utils.scale import scale_dem
import pyposeidon
import math


def get_hfun(dem, path=".", tag="jigsaw", resolution_min=0.01, resolution_max=0.5, dhdx=0.15, imax=100, **kwargs):

    # scale bathymetry
    try:
        b = dem.adjusted.to_dataframe()
    except:
        b = dem.elevation.to_dataframe()

    b = b.reset_index()
    b.columns = ["longitude", "latitude", "z"]

    nodes = scale_dem(b, resolution_min, resolution_max, **kwargs)

    x = dem.longitude.values
    y = dem.latitude.values

    tria = MakeTriangleFaces(y.shape[0], x.shape[0])

    points = np.column_stack([nodes.longitude, nodes.latitude])
    tria3 = pd.DataFrame(tria, columns=["a", "b", "c"])

    # Use Matplotlib for triangulation
    triang = matplotlib.tri.Triangulation(points[:, 0], points[:, 1], tria)
    #    tri3 = triang.triangles
    edges = triang.edges
    # edge lengths
    ptdiff = lambda p: (p[0][0] - p[1][0], p[0][1] - p[1][1])
    diffs = map(ptdiff, points[edges])
    elen = [math.hypot(d1, d2) for d1, d2 in diffs]

    hfun = nodes.d2.values
    hfun = hfun.reshape(hfun.shape[0], -1)

    [fun, flag] = limgrad2(edges, elen, hfun, dhdx, imax)

    cfun = fun.flatten().reshape((y.shape[0], x.shape[0]))

    dh = xr.Dataset(
        {"h": (["latitude", "longitude"], cfun)},
        coords={"longitude": ("longitude", x), "latitude": ("latitude", y)},
    )

    return dh


def to_global_hfun(nodes, elems, fpos, **kwargs):

    elems["d"] = 0
    nodes["z"] = 0

    R = kwargs.get("R", 1.0)

    sv = 4 * R ** 2 / (nodes.u ** 2 + nodes.v ** 2 + 4 * R ** 2)
    nodes["h"] = nodes.d2 / sv

    out = xr.merge([nodes.to_xarray(), elems.to_xarray()])

    to_hfun_mesh(out, fpos)

    ## make dataset
    els = xr.DataArray(
        elems.loc[:, ["a", "b", "c"]],
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    nod = (
        nodes.loc[:, ["u", "v"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "u": "SCHISM_hgrid_node_x",
                "v": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    bg = xr.Dataset({"h": (["nSCHISM_hgrid_node"], nodes.h.values)})

    dh = xr.merge([nod, els, bg])

    return dh


def to_hfun_mesh(dh, fhfun):

    dps = dh[["u", "v", "z"]].to_dataframe().dropna()
    hs = dh[["h"]].to_dataframe().dropna()
    trii = dh[["a", "b", "c", "d"]].to_dataframe().dropna()

    with open(fhfun, "w") as f:
        f.write("#{}; created by pyposeidon\n".format(pyposeidon.__version__))
        f.write("MSHID=3;EUCLIDEAN-MESH\n")
        f.write("NDIMS=2\n")
        f.write("POINT={}\n".format(dps.shape[0]))

    with open(fhfun, "a") as f:
        dps.to_csv(f, index=False, header=0, sep=";")

    with open(fhfun, "a") as f:
        f.write("VALUE={};1\n".format(dps.shape[0]))
        hs.to_csv(f, index=False, header=0)

    with open(fhfun, "a") as f:
        f.write("TRIA3={}\n".format(trii.shape[0]))
        trii.to_csv(f, index=False, header=0, sep=";")


def to_hfun_grid(dh, fhfun):
    # write hfun file

    # write header
    with open(fhfun, "w") as f:
        f.write("#{}; created by pyposeidon\n".format(pyposeidon.__version__))
        f.write("MSHID=3;EUCLIDEAN-GRID\n")
        f.write("NDIMS=2\n")
        f.write("COORD=1;{}\n".format(dh.longitude.size))

    with open(fhfun, "a") as f:
        np.savetxt(f, dh.longitude.values)

    with open(fhfun, "a") as f:
        f.write("COORD=2;{}\n".format(dh.latitude.size))

    with open(fhfun, "a") as f:
        np.savetxt(f, dh.latitude.values)

    with open(fhfun, "a") as f:
        f.write("VALUE={};1\n".format(dh.h.size))

    with open(fhfun, "a") as f:
        np.savetxt(f, dh.h.T.values.flatten())
