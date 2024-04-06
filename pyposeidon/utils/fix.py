"""
Mesh adjustment functions

"""

# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import geopandas as gp
import shapely
import pyresample
import pandas as pd
import xarray as xr
import sys
import os
import shutil
from glob import glob
from tqdm.auto import tqdm
from pyposeidon.utils.coastfix import simplify

# logging setup
import logging

logger = logging.getLogger(__name__)


def get_tiles(data, chunks=None):

    # chuck
    if chunks:
        ilats = data.elevation.chunk({"longitude": chunks[0], "latitude": chunks[1]}).chunks[0]
        ilons = data.elevation.chunk({"longitude": chunks[0], "latitude": chunks[1]}).chunks[1]
    else:
        ilats = data.elevation.chunk("auto").chunks[0]
        ilons = data.elevation.chunk("auto").chunks[1]

    if len(ilons) == 1:
        ilons = (int(ilons[0] / 2), int(ilons[0] / 2))

    idx = [sum(ilons[:i]) for i in range(len(ilons) + 1)]
    jdx = [sum(ilats[:i]) for i in range(len(ilats) + 1)]

    blon = list(zip(idx[:-1], idx[1:]))
    blat = list(zip(jdx[:-1], jdx[1:]))

    perms = [(x, y) for x in blon for y in blat]

    return perms


def fix_dem(dem, coastline, buffer=0.0, **kwargs):

    perms = get_tiles(dem)

    i = 0
    check = True

    if not os.path.exists("./fixtmp/"):
        os.makedirs("./fixtmp/")

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

        de_, check_, flag = fix(de, coastline, **kwargs)

        ide = de_.sel(latitude=slice(lat1, lat2)).sel(longitude=slice(lon1, lon2))

        ide.to_netcdf("./fixtmp/ide{:03d}.nc".format(i))
        i += 1

        check = check and check_

    ifiles = glob("./fixtmp/ide*")

    fdem = xr.open_mfdataset(ifiles)

    ##cleanup
    shutil.rmtree("./fixtmp")

    return fdem, check


def fix(dem, coastline, **kwargs):
    # ---------------------------------------------------------------------
    logger.info("adjust dem\n")
    # ---------------------------------------------------------------------

    ifunction = kwargs.get("resample_function", "nearest")
    reset_flag = kwargs.get("reset_flag", False)

    # define coastline
    try:
        shp = gp.GeoDataFrame.from_file(coastline)
    except:
        shp = gp.GeoDataFrame(coastline)

    sc = kwargs.get("simplify_coastlines", False)

    if sc:
        shp = simplify(shp)

    if "ival" in dem.data_vars:
        xp = dem.ilons.values
        yp = dem.ilats.values
    else:
        xp = dem.longitude.values
        yp = dem.latitude.values

    minlon = xp.min()
    maxlon = xp.max()
    minlat = yp.min()
    maxlat = yp.max()

    if xp.mean() < 0 and xp.min() < -180.0:
        flag = -1
    elif xp.mean() > 0 and xp.max() > 180.0:
        flag = 1
    else:
        flag = 0

    if flag == 1:
        block1 = shp.cx[minlon:180, minlat:maxlat].copy()
        block2 = shp.cx[-180 : (maxlon - 360.0), minlat:maxlat].copy()

        for idx, poly in block2.iterrows():
            block2.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x + 360.0, y), poly.geometry)

        block = pd.concat([block1, block2])

    elif flag == -1:
        block1 = shp.cx[minlon + 360 : 180, minlat:maxlat].copy()
        block2 = shp.cx[-180:maxlon, minlat:maxlat].copy()

        for idx, poly in block1.iterrows():
            block1.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x - 360.0, y), poly.geometry)

        block = pd.concat([block1, block2])

    else:
        block = shp.cx[minlon:maxlon, minlat:maxlat]

    try:
        block = gp.GeoDataFrame(geometry=list(block.unary_union.geoms))
    except:
        pass

    # ---------------------------------------------------------------------
    logger.debug("compute water and land\n")
    # ---------------------------------------------------------------------

    # create a polygon of the lat/lon window
    grp = shapely.geometry.Polygon([(minlon, minlat), (minlon, maxlat), (maxlon, maxlat), (maxlon, minlat)])

    grp = grp.buffer(0.5)  # buffer it to get also the boundary points

    try:
        g = block.unary_union.symmetric_difference(grp)  # get the diff
    except:
        g = grp  # no land

    t = gp.GeoDataFrame({"geometry": [g]}).explode(index_parts=True).droplevel(0)

    t["length"] = t["geometry"][:].length  # optional

    t = t.sort_values(by="length", ascending=0)  # use the length to list them
    t = t.reset_index(drop=True)

    t["in"] = gp.GeoDataFrame(geometry=[grp.buffer(0.001)] * t.shape[0]).contains(t)  # find the largest of boundaries

    try:
        idx = np.where(t["in"] == True)
        b = t.iloc[idx].geometry
    except:
        b = shapely.geometry.GeometryCollection()

    b = b.unary_union

    # define wet/dry
    water = b
    land = grp - b

    if (not land) | (not water):
        # ---------------------------------------------------------------------
        logger.debug("only water/land present...\n")
        # ---------------------------------------------------------------------

        if "ival" in dem.data_vars:
            dem = dem.assign(fval=dem.ival)

        else:
            dem = dem.assign(adjusted=dem.elevation)

        return dem, True, flag

    if "ival" in dem.data_vars:
        df = pd.DataFrame(
            {
                "longitude": dem.ilons.values.flatten(),
                "latitude": dem.ilats.values.flatten(),
                "elevation": dem.ival.values.flatten(),
            }
        )
    else:
        df = dem.elevation.to_dataframe().reset_index()

    # ---------------------------------------------------------------------
    logger.debug("invoke shapely\n")
    # ---------------------------------------------------------------------

    spoints_ = shapely.points(
        list(df.loc[:, ["longitude", "latitude"]].values)
    )  # create shapely objects for the points

    # Add land boundaries to a shapely object
    try:
        lbs = []
        for l in tqdm(range(len(land.boundary.geoms))):
            z = shapely.linearrings(land.boundary.geoms[l].coords[:])
            lbs.append(z)
    except:
        lbs = shapely.linearrings(land.boundary.coords[:])

    bp = shapely.polygons(lbs)

    # ---------------------------------------------------------------------
    logger.debug("find wet and dry masks\n")
    # ---------------------------------------------------------------------

    # find the points on land

    tree = shapely.STRtree(spoints_)

    try:
        wl = []
        for l in tqdm(range(len(land.boundary.geoms))):
            wl.append(tree.query(bp[l], predicate="contains").tolist())
        ns = [j for i in wl for j in i]
    except:
        wl = tree.query(bp, predicate="contains").tolist()
        ns = wl

    lmask = np.zeros(spoints_.shape, dtype=bool)
    lmask[ns] = True

    wmask = ~lmask  # invert for wet mask

    # ---------------------------------------------------------------------
    logger.debug("fix wet points\n")
    # ---------------------------------------------------------------------

    # Now see if the wet points have indeed negative values

    pw_mask = df.loc[wmask, "elevation"] > 0

    if pw_mask.sum() > 0:
        pw = df.loc[wmask][pw_mask]  # problematic points: bathymetry > 0 in wet area

        # Resample to fix that ...
        xw = pw.longitude.values
        yw = pw.latitude.values

        bw = resample(dem, xw, yw, var="elevation", wet=True, flag=flag, reset_flag=reset_flag, function=ifunction)

        df.loc[pw.index, "elevation"] = bw  # replace in original dataset

    # ---------------------------------------------------------------------
    logger.debug("fix dry points\n")
    # ---------------------------------------------------------------------

    # .. the same for dry points

    pl_mask = df.loc[lmask, "elevation"] < 0

    if pl_mask.sum() > 0:
        pl = df.loc[lmask][pl_mask]  # problematic points: bathymetry <0 in dry area

        ## Resample to fix that
        xl = pl.longitude.values
        yl = pl.latitude.values

        bd = resample(dem, xl, yl, var="elevation", wet=False, flag=flag, reset_flag=reset_flag, function=ifunction)

        df.loc[pl.index, "elevation"] = bd  # replace in original dataset

    # ---------------------------------------------------------------------
    logger.debug("assemble dataset \n")
    # ---------------------------------------------------------------------

    # reassemble dataset

    if "ival" in dem.data_vars:
        if len(dem.ival.shape) == 1:
            new_dem = df.elevation.to_xarray()
            new_dem = xr.merge([new_dem])
            new_dem = new_dem.rename({"elevation": "fval"})
            new_dem.fval.attrs = {"coastline": "based on coastline"}
            new_dem = new_dem.rename({"index": "k"}).drop_vars("k")
        else:
            new_dem = (
                df.set_index(["latitude", "longitude"])
                .to_xarray()
                .rename({"longitude": "l", "latitude": "k", "elevation": "fval"})
                .drop_vars(["k", "l"])
            )

    else:
        df_new = df.set_index(["latitude", "longitude"])
        new_dem = df_new.to_xarray()
        new_dem = new_dem.rename({"elevation": "adjusted"})
        new_dem.attrs = {"coastline": "based on coastline"}

    cdem = xr.merge([dem, new_dem])

    nanp = check1(cdem, water)

    if len(nanp) == 0:
        valid = True
    else:
        valid = False

    logger.info("Nan value for {} sea points".format(len(nanp)))

    check = kwargs.get("check", False)

    if check:
        on_coast = check2(cdem, shp)

        logger.info("{} points on the boundary, setting to zero".format(len(on_coast)))

        if "fval" in cdem.data_vars:
            tt = len(cdem.fval.shape)

            if tt == 1:
                cdem.fval[on_coast] = 0.0

            elif tt == 2:
                bmask = np.zeros(cdem.fval.shape, dtype=bool)  # create mask
                for idx, [i, j] in enumerate(on_coast):
                    bmask[i, j] = True
                cdem.fval.values[bmask] = 0.0  # set value

        elif "adjusted" in cdem.data_vars:
            bmask = np.zeros(cdem.adjusted.shape, dtype=bool)  # create mask
            for idx, [i, j] in enumerate(on_coast):
                bmask[i, j] = True
            cdem.adjusted.values[bmask] = 0.0  # set value

            logger.info(f"setting {bmask.size} land points with nan values to zero")
            cdem["adjusted"] = cdem.adjusted.fillna(0.0)  # for land points if any (lakes, etc.)

    else:
        logger.info("setting land points with nan values to zero")
        if "ival" in cdem.data_vars:
            cdem["fval"] = cdem.fval.fillna(0.0)
        elif "adjusted" in cdem.data_vars:
            cdem["adjusted"] = cdem.adjusted.fillna(0.0)

    return cdem, valid, flag


def check1(dataset, water):
    # check it

    if "fval" in dataset.data_vars:
        tt = len(dataset.fval.shape)

        if tt == 1:
            ids = np.argwhere(dataset.fval.values.flatten() > 0).flatten()

            xp = dataset.ilons.values.flatten()[ids]
            yp = dataset.ilats.values.flatten()[ids]

        elif tt == 2:
            xy = np.argwhere(dataset.fval.values > 0)

            xx = [x for [x, y] in xy]
            yy = [y for [x, y] in xy]

            smask = dataset.fval.values > 0

            xp = dataset.ilons.values[smask]
            yp = dataset.ilats.values[smask]

    elif "adjusted" in dataset.data_vars:
        xy = np.argwhere(dataset.adjusted.values > 0)

        coords = [0, 0]

        coords[0] = [x for [x, y] in xy]
        coords[1] = [y for [x, y] in xy]

        dims = list(dataset.adjusted.dims)

        match1 = dims.index("longitude")
        match2 = dims.index("latitude")

        xp = dataset.longitude[coords[match1]].values
        yp = dataset.latitude[coords[match2]].values

    wet = gp.GeoDataFrame(geometry=[water])

    nanpoints = [shapely.Point([x, y]) for (x, y) in zip(xp, yp)]

    tree = shapely.STRtree(nanpoints)

    wn = tree.query(water, predicate="contains").tolist()

    return wn


def check2(dataset, coastline):
    # check it

    tt = None

    if "fval" in dataset.data_vars:
        tt = len(dataset.fval.shape)

        if tt == 1:
            ids = np.argwhere(dataset.fval.values.flatten() > 0).flatten()

            xp = dataset.ilons.values.flatten()[ids]
            yp = dataset.ilats.values.flatten()[ids]

        elif tt == 2:
            xy = np.argwhere(dataset.fval.values > 0)

            xx = [x for [x, y] in xy]
            yy = [y for [x, y] in xy]

            smask = dataset.fval.values > 0

            xp = dataset.ilons.values[smask]
            yp = dataset.ilats.values[smask]

    elif "adjusted" in dataset.data_vars:
        tt = 0

        xy = np.argwhere(dataset.adjusted.values > 0)

        coords = [0, 0]

        coords[0] = [x for [x, y] in xy]
        coords[1] = [y for [x, y] in xy]

        dims = list(dataset.adjusted.dims)

        match1 = dims.index("longitude")
        match2 = dims.index("latitude")

        xp = dataset.longitude[coords[match1]].values
        yp = dataset.latitude[coords[match2]].values

    cpoints = [shapely.Point([x, y]) for (x, y) in zip(xp, yp)]

    xps = gp.GeoDataFrame(geometry=cpoints)

    tree = shapely.STRtree(xps.buffer(0.0001).geometry.values)

    if not coastline.boundary.is_empty.all():
        coasts = gp.GeoDataFrame(geometry=coastline.boundary)
    else:
        coasts = coastline

    cc = coasts.unary_union

    wn = tree.query(cc, predicate="intersects").tolist()

    wn.sort()

    if tt == 0:
        bps = xy[wn].tolist()

    elif tt == 1:
        bps = ids[wn]

    elif tt == 2:
        bps = [[xx[i], yy[i]] for i in wn]

    elif tt == None:
        bps = None

    return bps


def resample(dem, xw, yw, var=None, wet=True, flag=None, reset_flag=False, function="nearest"):
    # Define points with positive bathymetry
    x, y = np.meshgrid(dem.longitude, dem.latitude)

    print(f"reset_flag={reset_flag}")

    if reset_flag:
        flag = 0

    if flag == 1:
        gx = xw - 180.0
        xx = x - 180.0
    elif flag == -1:
        gx = xw + 180.0
        xx = x + 180.0
    else:
        gx = xw
        xx = x

    # fill the nan, if present, with values in order to compute values there if needed.
    dem[var].data[np.isnan(dem[var].values)] = 9999.0

    if wet:
        mx = np.ma.masked_array(xx, dem[var].values > 0)
        my = np.ma.masked_array(y, dem[var].values > 0)

        # mask positive bathymetry
        mdem = np.ma.masked_array(dem[var], dem[var].values > 0)

    else:
        mx = np.ma.masked_array(xx, dem[var].values < 0)
        my = np.ma.masked_array(y, dem[var].values < 0)

        # mask positive bathymetry
        mdem = np.ma.masked_array(dem[var], dem[var].values < 0)

    orig = pyresample.geometry.SwathDefinition(lons=mx, lats=my)  # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=gx, lats=yw)  # wet points

    mdem = mdem.astype(float)

    if function == "nearest":
        bw = pyresample.kd_tree.resample_nearest(orig, mdem, targ, radius_of_influence=100000, fill_value=np.nan)

    elif function == "gauss":
        bw = pyresample.kd_tree.resample_gauss(
            orig, mdem, targ, radius_of_influence=500000, neighbours=10, sigmas=250000, fill_value=np.nan
        )

    return bw


def dem_range(data, lon_min, lon_max, lat_min, lat_max):
    dlon0 = round(data.longitude.data.min())
    dlon1 = round(data.longitude.data.max())

    # recenter the window
    if dlon1 - dlon0 == 360.0:
        lon0 = lon_min + 360.0 if lon_min < data.longitude.min() else lon_min
        lon1 = lon_max + 360.0 if lon_max < data.longitude.min() else lon_max

        lon0 = lon0 - 360.0 if lon0 > data.longitude.max() else lon0
        lon1 = lon1 - 360.0 if lon1 > data.longitude.max() else lon1

    else:
        lon0 = lon_min
        lon1 = lon_max

    if (lon_min < data.longitude.min()) or (lon_max > data.longitude.max()):
        logger.info("Lon must be within {} and {}".format(data.longitude.min().values, data.longitude.max().values))
        logger.info("compensating if global dataset available")

    if (lat_min < data.latitude.min()) or (lat_max > data.latitude.max()):
        logger.info("Lat is within {} and {}".format(data.latitude.min().values, data.latitude.max().values))

    # get idx
    if lon_max - lon_min == dlon1 - dlon0:
        i0 = 0 if lon_min == dlon0 else int(data.longitude.shape[0] / 2) + 2  # compensate for below
        i1 = data.longitude.shape[0] if lon_max == dlon1 else -int(data.longitude.shape[0] / 2) - 2
    else:
        i0 = np.abs(data.longitude.data - lon0).argmin()
        i1 = np.abs(data.longitude.data - lon1).argmin()

    j0 = np.abs(data.latitude.data - lat_min).argmin()
    j1 = np.abs(data.latitude.data - lat_max).argmin()

    # expand the window a little bit
    lon_0 = max(0, i0 - 2)
    lon_1 = min(data.longitude.size, i1 + 2)

    lat_0 = max(0, j0 - 2)
    lat_1 = min(data.latitude.size, j1 + 2)

    # descenting lats
    if j0 > j1:
        j0, j1 = j1, j0
        lat_0 = max(0, j0 - 1)
        lat_1 = min(data.latitude.size, j1 + 3)

    if i0 > i1:
        p1 = data.isel(longitude=slice(lon_0, data.longitude.size), latitude=slice(lat_0, lat_1))

        p1 = p1.assign_coords({"longitude": p1.longitude.values - 360.0})

        p2 = data.isel(longitude=slice(0, lon_1), latitude=slice(lat_0, lat_1))

        dem = xr.concat([p1, p2], dim="longitude")

    else:
        dem = data.isel(longitude=slice(lon_0, lon_1), latitude=slice(lat_0, lat_1))

    if np.abs(np.mean(dem.longitude) - np.mean([lon_min, lon_max])) > 170.0:
        c = np.sign(np.mean([lon_min, lon_max]))
        dem["longitude"] = dem["longitude"] + c * 360.0

    dem_data = xr.merge([dem])

    return dem_data
