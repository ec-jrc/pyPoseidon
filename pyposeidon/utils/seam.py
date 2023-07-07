import pandas as pd
import geopandas as gp
import numpy as np
import shapely
from glob import glob
import os
from tqdm.auto import tqdm
import pyresample
import xarray as xr
from shapely.ops import triangulate
import logging

logger = logging.getLogger(__name__)


def get_seam(x, y, z, tri3, **kwargs):
    if z == None:
        z = 1.0

    # store lon/lat/elems
    gr = pd.DataFrame({"lon": x, "lat": y, "z": z})
    elems = pd.DataFrame(tri3, columns=["a", "b", "c"])

    # discover crossover elements
    a, b, c = tri3.T
    lon = x
    lon_a = lon[a]
    lon_b = lon[b]
    lon_c = lon[c]

    max_lon = kwargs.get("max_lon", 340)
    idl = np.where(
        ((lon_a * lon_b < 0) & (np.abs(lon_a - lon_b) >= max_lon))
        | ((lon_a * lon_c < 0) & (np.abs(lon_a - lon_c) >= max_lon))
        | ((lon_b * lon_c < 0) & (np.abs(lon_b - lon_c) >= max_lon))
    )[0]

    mels = elems.iloc[idl].copy().reset_index(drop=True)

    # polygonize
    mels["ap"] = mels.apply(lambda x: gr.loc[x.a, ["lon", "lat"]].values, axis=1)
    mels["bp"] = mels.apply(lambda x: gr.loc[x.b, ["lon", "lat"]].values, axis=1)
    mels["cp"] = mels.apply(lambda x: gr.loc[x.c, ["lon", "lat"]].values, axis=1)

    mels["geometry"] = mels.apply(lambda x: shapely.geometry.Polygon([x.ap, x.bp, x.cp]), axis=1)

    # Split crossing elements
    adv = 0
    p1 = []
    p2 = []
    for indx, vals in mels.iterrows():
        #    print(indx)
        pol = vals.geometry
        lons = [x for (x, y) in list(pol.exterior.coords)]
        incl = np.array(lons[:-1]).mean()
        if incl < 0.0:
            lons = [x - 360.0 if x > 0.0 else x for x in pol.boundary.xy[0]]
        else:
            lons = [x + 360.0 if x < 0.0 else x for x in pol.boundary.xy[0]]

        nns = list(zip(lons, list(pol.boundary.xy[1])))  # points of recasted element

        npol = shapely.geometry.LineString(nns)  # make a geometrical object (triangle) out of the nodes above

        # create a meridian line
        if incl < 0.0:
            l = shapely.geometry.LineString([[-180.0, -90.0], [-180.0, 90.0]])
        else:
            l = shapely.geometry.LineString([[180.0, -90.0], [180.0, 90.0]])

        cp = npol.intersection(l)  # find the intersection with the element

        try:
            cpp = [(x.coords[0][0], x.coords[0][1]) for x in cp.geoms]  # get cross nodes
        except:
            cpp = []

        de = pd.DataFrame(cpp + nns[:-1], columns=["lon", "lat"])

        de = de.drop_duplicates().reset_index(drop=True)

        points = shapely.geometry.MultiPoint(nns[:-1] + cpp)

        triangles = triangulate(points)

        tes = []
        for triangle in triangles:
            blon = [x for (x, y) in list(triangle.exterior.coords)]
            blat = [y for (x, y) in list(triangle.exterior.coords)]
            tes.append(de[(de["lon"].isin(blon)) & (de["lat"].isin(blat))].index.values)

        nels = pd.DataFrame(tes, columns=["a", "b", "c"])

        if cpp:
            de = pd.concat([de, de.loc[:1]], ignore_index=True)  # replicate the meridian cross points
            de.lon[-2:] *= -1

            if de[de.lon > 180.0].size > 0:
                ids = de[de.lon > 180.0].index.values  # problematic nodes
                de.loc[de.lon > 180.0, "lon"] -= 360.0
            elif de[de.lon < -180.0].size > 0:
                ids = de[de.lon < -180.0].index.values  # problematic nodes
                de.loc[de.lon < -180.0, "lon"] += 360.0

        p1.append(de)

        if cpp:
            des = nels.loc[(nels.a.isin(ids)) | (nels.b.isin(ids)) | (nels.c.isin(ids))].copy()
            des.loc[des.a == 0, "a"] = de.index[-2]
            des.loc[des.a == 1, "a"] = de.index[-1]
            des.loc[des.b == 0, "b"] = de.index[-2]
            des.loc[des.b == 1, "b"] = de.index[-1]
            des.loc[des.c == 0, "c"] = de.index[-2]
            des.loc[des.c == 1, "c"] = de.index[-1]
            nels.loc[des.index] = des

        nels = nels + adv  # reindex to global index

        p2.append(nels)

        adv = nels.values.max() + 1

    ng = pd.concat(p1)
    ng["z"] = 1  ## ADJUST

    ng.reset_index(inplace=True, drop=True)

    nge = pd.concat(p2)
    nge.reset_index(inplace=True, drop=True)

    si = gr.index[-1]
    ## drop the problematic elements
    ges = elems.drop(idl)
    ## append new nodes
    mes = pd.concat([gr, ng])

    ## Make new elements index global
    nges = nge + si + 1
    # append new elements
    ges = pd.concat([ges, nges])
    ges.reset_index(inplace=True, drop=True)
    mes.reset_index(inplace=True, drop=True)

    xx = mes.lon.values
    yy = mes.lat.values
    return xx, yy, ges.values


def to_2d(dataset=None, var=None, mesh=None, **kwargs):
    x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
    y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
    tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")

    x = dataset[x_var][:].values
    y = dataset[y_var][:].values
    tri3 = dataset[tes_var].values[:, :3].astype(int)

    if mesh is not None:
        [xn, yn, tri3n] = mesh
    else:
        xn, yn, tri3n = get_seam(x, y, None, tri3)

    nps = xn.shape[0] - x.shape[0]  # number of extra nodes
    xi = xn[-nps:].copy()  # lat/lon of extra nodes
    yi = yn[-nps:].copy()

    # reposition to avoid -180/180 boundary
    pxi = reposition(xi)
    pyi = yi

    a = pxi.min()
    b = pxi.max() - 360

    xmask = (x < b - 10) | (x > a + 10)  # get a slice of the mesh around IM
    # reposition to avoid -180/180 boundary
    px = reposition(x[xmask])
    py = y[xmask]

    d = np.mean([px.max() - px.min(), pxi.max() - pxi.min()])

    # create swaths
    orig = pyresample.geometry.SwathDefinition(lons=px - d, lats=py)
    targ = pyresample.geometry.SwathDefinition(lons=pxi - d, lats=pyi)

    # Resample
    if len(dataset[var].shape) == 1:
        z = dataset[var].values
        zm = z[xmask]
        z_ = pyresample.kd_tree.resample_nearest(orig, zm, targ, radius_of_influence=200000)  # , fill_value=0)
        xelev = np.concatenate((z, z_))

        # create xarray
        xe = xr.Dataset(
            {
                var: (["nSCHISM_hgrid_node"], xelev),
                "SCHISM_hgrid_node_x": (["nSCHISM_hgrid_node"], xn),
                "SCHISM_hgrid_node_y": (["nSCHISM_hgrid_node"], yn),
                "SCHISM_hgrid_face_nodes": (
                    ["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
                    tri3n,
                ),
            },
        )

    elif "time" in dataset[var].coords:
        it_start = kwargs.get("it_start", 0)
        it_end = kwargs.get("it_end", dataset.time.shape[0])

        if not os.path.exists("./seamtmp/"):
            os.makedirs("./seamtmp/")

        for i in tqdm(range(it_start, it_end)):
            z = dataset[var].values[i, :]
            zm = z[xmask]
            z_ = pyresample.kd_tree.resample_nearest(orig, zm, targ, radius_of_influence=200000, fill_value=0)
            e = np.concatenate((z, z_))
            e = e[np.newaxis, :]  # make 2d

            # create xarray
            xi = xr.Dataset(
                {
                    var: (["time", "nSCHISM_hgrid_node"], e),
                    "SCHISM_hgrid_node_x": (["nSCHISM_hgrid_node"], xn),
                    "SCHISM_hgrid_node_y": (["nSCHISM_hgrid_node"], yn),
                    "SCHISM_hgrid_face_nodes": (
                        ["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
                        tri3n,
                    ),
                },
                coords={"time": ("time", [dataset.time.values[i]])},
            )

            xi.to_netcdf("./seamtmp/x_{:03d}.nc".format(i))

        xe = xr.open_mfdataset("./seamtmp/x_*.nc", data_vars="minimal")

        # cleanup
        xfiles = glob("./seamtmp/x_*.nc")
        for f in xfiles:
            os.remove(f)
        os.removedirs("./seamtmp/")

    return xe


def reposition(px):
    px[px < 0] = px[px < 0] + 360.0

    return px
