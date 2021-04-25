import pandas as pd
import numpy as np
import geopandas as gp
import logging
import shapely
import sys

logger = logging.getLogger("pyposeidon")


def tag_(**kwargs):

    world = kwargs.get("coastlines", None)

    if world is None:
        logger.error("coastlines not given")
        sys.exit(1)

    world = world.explode()

    geometry = kwargs.get("geometry", None)

    try:
        lon_min = geometry["lon_min"]
        lon_max = geometry["lon_max"]
        lat_min = geometry["lat_min"]
        lat_max = geometry["lat_max"]
    except:
        logger.error("geometry not set properly")
        sys.exit(1)

    # create a polygon of the lat/lon window
    grp = shapely.geometry.Polygon([(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_max), (lon_max, lat_min)])

    # create a LineString of the grid
    grl = shapely.geometry.LineString(
        [(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_max), (lon_max, lat_min), (lon_min, lat_min)]
    )

    # check -180/180 trespass
    if np.mean([lon_min, lon_max]) < 0 and lon_min < -180.0:
        flag = -1
    elif np.mean([lon_min, lon_max]) > 0 and lon_max > 180.0:
        flag = 1
    else:
        flag = 0

    # adjust abd mask based on lat/lon window
    if flag == 1:
        block1 = world.cx[lon_min:180, lat_min:lat_max].copy()
        block2 = world.cx[-180 : (lon_max - 360.0), lat_min:lat_max].copy()

        for idx, poly in block2.iterrows():
            block2.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x + 360.0, y), poly.geometry)

        block = block1.append(block2)

    elif flag == -1:

        block1 = world.cx[lon_min + 360 : 180, lat_min:lat_max].copy()
        block2 = world.cx[-180:lon_max, lat_min:lat_max].copy()

        for idx, poly in block1.iterrows():
            block1.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x - 360.0, y), poly.geometry)

        block = block1.append(block2)

    else:
        block = world.cx[lon_min:lon_max, lat_min:lat_max]

    #    block = block.to_crs('epsg=3763').buffer(10).to_crs("EPSG:4326")

    g = block.unary_union.symmetric_difference(grp)  # get the dif from the world

    try:  # make geoDataFrame
        t = gp.GeoDataFrame({"geometry": g})
    except:
        t = gp.GeoDataFrame({"geometry": [g]})
    t["length"] = t["geometry"][:].length  # get length
    t = t.sort_values(by="length", ascending=0)  # sort
    t = t.reset_index(drop=True)

    t["in"] = gp.GeoDataFrame(geometry=[grp] * t.shape[0]).contains(t)  # find the largest of boundaries
    idx = np.where(t["in"] == True)[0][0]  # first(largest) boundary within lat/lon

    b = t.iloc[idx].geometry  # get the largest

    # SETUP JIGSAW

    dic = {}
    try:
        for l in range(len(b.boundary)):
            lon = []
            lat = []
            for x, y in b.boundary[l].coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({"line{}".format(l): {"lon": lon, "lat": lat}})
    except:
        lon = []
        lat = []
        for x, y in b.boundary.coords[:]:
            lon.append(x)
            lat.append(y)
        dic.update({"line{}".format(0): {"lon": lon, "lat": lat}})

    dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}
    df = pd.concat(dict_of_df, axis=0)
    df["z"] = 0
    df = df.drop_duplicates()  # drop the repeat value on closed boundaries

    # open (water) boundaries

    try:
        water = b.boundary[0] - (b.boundary[0] - grl)
    except:
        water = b.boundary - (b.boundary - grl)

    try:
        cwater = shapely.ops.linemerge(water)
    except:
        cwater = water

    mindx = 1
    # get all lines in a pandas DataFrame
    if cwater.type == "LineString":
        lon = []
        lat = []
        for x, y in cwater.coords[:]:
            lon.append(x)
            lat.append(y)
        dic = {"line{}".format(mindx): {"lon": lon, "lat": lat, "z": 0, "tag": mindx}}
        mindx += 1

    elif cwater.type == "MultiLineString":
        dic = {}
        for l in range(len(cwater)):
            lon = []
            lat = []
            for x, y in cwater[l].coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({"line{}".format(mindx): {"lon": lon, "lat": lat, "z": 0, "tag": mindx}})
            mindx += 1

    dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}
    df_water = pd.concat(dict_of_df, axis=0)

    #    ibs = df_water.index.levels[0].shape[0]
    #    if ibs > 0 :

    ## resample boundaries
    #        ndfs={}
    #        for ic in tqdm(range(ibs)):
    #            contour=df_water.index.levels[0][ic]
    #            curve=df_water.loc[contour,['lon','lat']]
    #            if curve.shape[0] < 100:
    #                di = spline(curve, 100)
    #                di['z']=df_water.loc[contour].z.values[0]
    #                di['tag']=df_water.loc[contour].tag.values[0].astype(int)
    #            else:
    #                di = curve
    #            ndfs.update({contour:di})

    #        df_water = pd.concat(ndfs, axis=0)

    #        df_water['z'] = df_water.z.values.astype(int)
    #        df_water['tag'] = df_water.tag.values.astype(int)

    # land boundaries!!
    try:
        land = b.boundary[0] - grl
    except:
        land = b.boundary - grl

    try:
        cland = shapely.ops.linemerge(land)
    except:
        cland = land

    mindx = -1

    # get all lines in a pandas DataFrame
    dic_land = {}

    if cland.type == "LineString":
        lon = []
        lat = []
        for x, y in cland.coords[:]:
            lon.append(x)
            lat.append(y)
        dic_land = {"line{}".format(mindx): {"lon": lon, "lat": lat, "z": 0, "tag": mindx}}
        mindx -= 1

    elif cland.type == "MultiLineString":
        dic_land = {}
        for l in range(len(cland)):
            lon = []
            lat = []
            for x, y in cland[l].coords[:]:
                lon.append(x)
                lat.append(y)
            dic_land.update({"line{}".format(mindx): {"lon": lon, "lat": lat, "z": 0, "tag": mindx}})
            mindx -= 1

    dict_of_df = {k: pd.DataFrame(v) for k, v in dic_land.items()}

    try:
        df_land = pd.concat(dict_of_df, axis=0)
    except:
        df_land = pd.DataFrame({})

    ### interpolate if needed

    even = kwargs.get("even", False)
    ds = kwargs.get("ds", 0.001)

    if even:

        ibs = df_land.index.levels[0].shape[0]

        if ibs < 0:
            ## resample boundaries
            ndfs = {}
            for ic in tqdm(range(ibs)):
                contour = df_land.index.levels[0][ic]
                curve = df_land.loc[contour, ["lon", "lat"]]
                di = spline(curve, ds=ds)
                di["z"] = df_land.loc[contour].z.values[0]
                di["tag"] = df_land.loc[contour].tag.values[0].astype(int)
                ndfs.update({contour: di})

            df_land = pd.concat(ndfs, axis=0)
            df_land["z"] = df_land.z.values.astype(int)
            df_land["tag"] = df_land.tag.values.astype(int)

    # put together water & land
    ddf = pd.concat([df_water, df_land])

    # Sort outer boundary

    out_b = []
    for line in ddf.index.levels[0]:
        out_b.append(shapely.geometry.LineString(ddf.loc[line, ["lon", "lat"]].values))

    merged = shapely.ops.linemerge(out_b)
    merged = pd.DataFrame(merged.coords[:], columns=["lon", "lat"])
    merged = merged.drop_duplicates()
    match = ddf.drop_duplicates(["lon", "lat"]).droplevel(0)
    match = match.reset_index(drop=True)

    df1 = merged.sort_values(["lon", "lat"])
    df2 = match.sort_values(["lon", "lat"])
    df2.index = df1.index
    final = df2.sort_index()
    final = pd.concat([final], keys=["line0"])

    bmindx = mindx

    # merge with islands
    ndf = df.drop("line0")
    ndf["tag"] = ""
    for line in ndf.index.levels[0][1:]:
        ndf.loc[line, "tag"] = mindx - 1
        mindx -= 1
    ndf["tag"] = ndf.tag.astype(int)

    ### interpolate if needed
    if even:

        ibs = ndf.index.levels[0].shape[0]

        ## resample boundaries
        ndfs = {}
        for ic in tqdm(range(ibs)):
            contour = ndf.index.levels[0][ic]
            curve = ndf.loc[contour, ["lon", "lat"]]
            curve = pd.concat([curve, curve.loc[0:0]]).reset_index(
                drop=True
            )  # add the first point to do a better interpolation
            di = spline(curve, ds=ds)
            di["z"] = ndf.loc[contour].z.values[0]
            di["tag"] = ndf.loc[contour].tag.values[0].astype(int)
            ndfs.update({contour: di.drop_duplicates(["lon", "lat"])})  # remove duplicated points (see above)

        ndf = pd.concat(ndfs, axis=0)

        ndf["z"] = ndf.z.values.astype(int)
        ndf["tag"] = ndf.tag.values.astype(int)

    df = pd.concat([final, ndf])

    return df, bmindx
