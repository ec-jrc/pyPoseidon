"""
Geometry module of pyposeidon. It manages model boundaries.
"""
from pyposeidon.utils.stereo import to_stereo
import pyposeidon.utils.pplot
import pyposeidon.dem as pdem
import pandas as pd
import numpy as np
import geopandas as gp
import logging
import shapely
from tqdm import tqdm
import sys
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


def simplify(geo):

    logger.info("simplify coastalines dataset")
    geo = geo.loc[:, ["geometry"]].explode(index_parts=True).droplevel(0).reset_index(drop=True)

    # Just in case merge
    if (geo.geom_type == "LineString").all():
        geo_ = gp.GeoDataFrame(geometry=[shapely.ops.linemerge(geo.geometry.values)])
        geo = geo_.explode(index_parts=True).droplevel(0).reset_index(drop=True)

    if (geo.geom_type == "Polygon").all():
        try:
            geo_ = list(geo.buffer(0).unary_union)
        except TypeError:
            geo_ = [geo.buffer(0).unary_union]

        geo = gp.GeoDataFrame(geometry=geo_)
        geo = gp.GeoDataFrame(geometry=geo.buffer(0))

    return geo


class Boundary:
    def __init__(self, **kwargs):
        """
        Set model boundaries

        !!! danger ""
            Due to a limitation of the Library rendering the docstrings, all arguments are marked
            as `required`, nevertheless they are all `Optional` except geometry.

        Args:
            geometry (Union[dict, str, GeoDataFrame]): A `GeoDataFrame` or the path to a shapefile or
                a dict defining the lat/lon window.
            coastlines (Union[str, GeoDataFrame]): A `GeoDataFrame` or the path to a shapefile which
                describes the coastlines. Defaults to `None`.
            cbuffer (float): The buffer in arcs for extending the coastlines. Defaults to `None`.
            levels (list[floats]): The range of DEM values for extracting the boundaries.
                When one valus is present it defines inner coastlines. When two values exist they define
                the extent. Defaults to `None`.
            dem_source (str): Path or url to bathymetric data.
        """

        geometry = kwargs.get("geometry", None)
        coastlines = kwargs.get("coastlines", None)
        cbuffer = kwargs.get("cbuffer", None)
        blevels = kwargs.get("blevels", None)
        prad = kwargs.get("R", 1.0)

        # COASTLINES
        if coastlines is None:
            logger.warning("coastlines not given")
            self.coasts = None

        elif isinstance(coastlines, str):
            logger.info("reading {}".format(coastlines))
            coasts = gp.GeoDataFrame.from_file(coastlines)
            self.coasts = simplify(coasts)
        elif isinstance(coastlines, gp.GeoDataFrame):
            logger.warning("coastlines is not a file, trying with geopandas Dataset")
            try:
                self.coasts = simplify(coastlines)
            except:
                logger.error("coastlines argument not valid ")
                sys.exit(1)

        # GEOMETRY
        if geometry is None:
            if levels is None:
                logger.error("geometry nor levels is given, exiting")
                sys.exit(1)

        if isinstance(geometry, dict):

            if self.coasts is None:
                logger.warning("coastlines might be required")

            self.geometry = geometry

        elif isinstance(geometry, str):

            if geometry == "global":

                if self.coasts is None:
                    logger.warning("coastlines might be required")

                self.geometry = "global"

            else:

                try:
                    self.geometry = gp.read_file(geometry)
                except:
                    logger.warning("geometry is not a file, trying with geopandas Dataset")
                    if isinstance(geometry, gp.GeoDataFrame):
                        self.geometry = geometry
                    else:
                        logger.error("geometry argument not valid ")
                        sys.exit(1)

        else:

            try:
                self.geometry = gp.read_file(geometry)
            except:
                logger.warning("geometry is not a file, trying with geopandas Dataset")
                if isinstance(geometry, gp.GeoDataFrame):
                    self.geometry = geometry
                else:
                    logger.error("geometry argument not valid ")
                    sys.exit(1)

        # Define internal boundary as isovalue of DEM
        if blevels:
            dsource = kwargs.get("dem_source", None)
            if dsource is None:
                logger.error("dem_source is required")

            dem = pdem.Dem(geometry=self.geometry, dem_source=dsource)
            dem_ = dem.Dataset

            self.coasts = get_dem_contours(blevels, dem_)

        # get boundaries
        if isinstance(self.geometry, dict):
            df = tag(self.geometry, self.coasts, cbuffer, blevels)
        elif isinstance(self.geometry, str):
            if self.coasts is None:
                logger.error("coastlines are missing .. exiting\n")
                sys.exit(1)

            df = global_tag(self.coasts, cbuffer, blevels, R=prad)
        elif isinstance(self.geometry, gp.GeoDataFrame):
            df = self.geometry

        # line tag
        df.loc[df.tag == "island", "lindex"] = np.arange(-df[df.tag == "island"].shape[0], 0)
        df.loc[df.tag == "land", "lindex"] = 1000 + np.arange(1, df[df.tag == "land"].shape[0] + 1)
        df.loc[df.tag == "open", "lindex"] = np.arange(1, df[df.tag == "open"].shape[0] + 1)
        df = df.sort_values("lindex", ascending=False)
        df.lindex = df.lindex.astype(int)

        # number of points
        df["nps"] = df.apply(lambda row: len(row.geometry.xy[1]) - 1, axis=1)

        self.contours = df.reset_index(drop=True)

    def show(self):

        return self.contours.plot(
            column="tag",
            legend=True,
            legend_kwds={
                "loc": "upper center",
                "bbox_to_anchor": (0.5, 1.15),
                "ncol": 3,
                "fancybox": True,
                "shadow": True,
            },
        )


def buffer_(coasts, cbuffer):

    # check
    if coasts.empty:
        return coasts

    # shrink the world so not to get over
    limit = shapely.geometry.Polygon([[-180, -90], [180, -90], [180, 90], [-180, 90], [-180, -90]]).buffer(-cbuffer)
    de = coasts.intersection(limit).explode(index_parts=True).droplevel(0)

    # buffer
    w_ = gp.GeoDataFrame(geometry=de.buffer(cbuffer))  # in arcs

    # drop empty objects
    empty = w_.is_empty
    w_ = w_.loc[~empty]
    w_ = w_.reset_index(drop=True)

    # evaluate multiple boundaries
    mls = []
    for pos, pol in w_.itertuples():
        bb = pol.boundary
        try:
            len(bb)
            mls.append(pos)
        except:
            pass

    # keep only exterior
    for idx, val in w_.loc[mls].iterrows():
        b = shapely.geometry.Polygon(val.geometry.exterior)
        w_.loc[idx] = b

    # join
    wu = w_.unary_union
    wu = gp.GeoDataFrame(geometry=[wu]).explode(index_parts=True).droplevel(0).reset_index(drop=True)

    rings = wu.boundary.is_ring  # get multiple boundaries instance
    wi = wu.loc[~rings]

    # get exterior (if any)
    for idx, val in wi.iterrows():
        b = shapely.geometry.Polygon(val.geometry.exterior)
        wu.loc[idx] = b

    return wu


def get_dem_contours(blevels, dem):

    if len(blevels) == 1:
        dvalues = [blevels[0] - 1, blevels[0]]
    elif len(blevels) == 2:
        dvalues = blevels
    else:
        logger.error("blevels not set properly")
        sys.exit(1)

    # check if in range
    r_check = [dem.elevation.values.min() < x < dem.elevation.values.max() for x in dvalues]

    if not any(r_check):
        logger.warning(
            "blevels not within DEM range {}:{}".format(dem.elevation.values.min(), dem.elevation.values.max())
        )
        return gp.GeoDataFrame(geometry=[])
    try:
        dem.adjusted.plot.contour(levels=dvalues)
    except:
        dem.elevation.plot.contour(levels=dvalues)

    # get graph
    cs = plt.gca()

    # outer contour
    outc = []
    outl = cs.collections[0].get_paths()
    for l in range(len(outl)):
        v = outl[l].vertices
        try:
            outc.append(shapely.geometry.LineString(v))
        except:
            pass

    # inner contour
    inc = []
    inl = cs.collections[1].get_paths()
    for l in range(len(inl)):
        v = inl[l].vertices
        try:
            inc.append(shapely.geometry.LineString(v))
        except:
            pass

    # Get continuous area
    lon_min = dem.longitude.min()
    lon_max = dem.longitude.max()
    lat_min = dem.latitude.min()
    lat_max = dem.latitude.max()

    # create a polygon of the lat/lon window
    grp = shapely.geometry.Polygon([(lon_min, lat_min), (lon_min, lat_max), (lon_max, lat_max), (lon_max, lat_min)])

    # outer contour
    b1 = []
    for geo_ in outc:
        ee = shapely.ops.split(grp, geo_)
        try:
            b1.append(ee[1])
        except:
            pass
    # inner contour
    b2 = []
    for geo_ in inc:
        ee = shapely.ops.split(grp, geo_)
        try:
            b2.append(ee[1])
        except:
            pass

    w1 = gp.GeoDataFrame(geometry=b1)
    w1_ = w1.unary_union
    w1 = gp.GeoDataFrame(geometry=[w1_])

    w2 = gp.GeoDataFrame(geometry=b2)
    w2_ = w2.unary_union
    w2 = gp.GeoDataFrame(geometry=[w2_])

    plt.close("all")
    return w1.loc[~w1.is_empty]  # cleanup


def tag(geometry, coasts, cbuffer, blevels):

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
        [
            (lon_min, lat_min),
            (lon_min, lat_max),
            (lon_max, lat_max),
            (lon_max, lat_min),
            (lon_min, lat_min),
        ]
    )

    # check -180/180 trespass
    if np.mean([lon_min, lon_max]) < 0 and lon_min < -180.0:
        flag = -1
    elif np.mean([lon_min, lon_max]) > 0 and lon_max > 180.0:
        flag = 1
    else:
        flag = 0

    try:

        # adjust and mask based on lat/lon window
        if flag == 1:
            block1 = coasts.cx[lon_min:180, lat_min:lat_max].copy()
            block2 = coasts.cx[-180 : (lon_max - 360.0), lat_min:lat_max].copy()

            for idx, poly in block2.iterrows():
                block2.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x + 360.0, y), poly.geometry)

            block = pd.concat([block1, block2])

        elif flag == -1:

            block1 = coasts.cx[lon_min + 360 : 180, lat_min:lat_max].copy()
            block2 = coasts.cx[-180:lon_max, lat_min:lat_max].copy()

            for idx, poly in block1.iterrows():
                block1.loc[idx, "geometry"] = shapely.ops.transform(lambda x, y, z=None: (x - 360.0, y), poly.geometry)

            block = pd.concat([block1, block2])

        else:
            block = coasts.cx[lon_min:lon_max, lat_min:lat_max]

        # Fix polygons around international line
        if flag != 0:

            wc = block[(block.bounds.maxx == flag * 180) | (block.bounds.minx == flag * 180)]

            cs = []  # adjust values around zero (in projection - international line in Platee Carree)
            for idx, line in wc.itertuples():
                try:
                    x_ = [
                        flag * 180 if abs(x - flag * 180) < 1.0e-3 else x for (x, y) in line.exterior.coords[:]
                    ]  # polygons
                    y_ = [y for (x, y) in line.exterior.coords[:]]
                except:
                    x_ = [
                        flag * 180 if abs(x - flag * 180) < 1.0e-3 else x for (x, y) in line.coords[:]
                    ]  # LineStrings
                    y_ = [y for (x, y) in line.coords[:]]

                cs.append(shapely.geometry.Polygon(list(zip(x_, y_))))

            ww = gp.GeoDataFrame(geometry=cs)

            try:
                gw = gp.GeoDataFrame(
                    geometry=list(ww.buffer(0).unary_union)
                )  # merge the polygons that are split (around -180/180)
            except:
                gw = gp.GeoDataFrame(geometry=list(ww.values))

            if wc.geom_type.all() != "Polygon":
                gw = gp.GeoDataFrame(geometry=gw.boundary.values)

            block = pd.concat([block[~block.index.isin(wc.index)], gw])

    except:
        if coasts == None:
            block = gp.GeoDataFrame(geometry=[])
            pass

    # buffer the coastlines on demand
    if cbuffer:
        logger.info("Buffering...")
        block = buffer_(block, cbuffer)
        logger.info("...done")

    # polygonize if need be for getting the symetric difference
    if (block.geom_type == "LineString").all():
        gg = shapely.ops.polygonize_full(block.geometry.values)
        block = gp.GeoDataFrame(geometry=list(gg[0])).explode(index_parts=True).droplevel(0).reset_index(drop=True)

    # bypass blocks in case of isodem
    if blevels:
        block = coasts.copy()

    if not block.empty:
        g = block.unary_union.symmetric_difference(grp)  # get the dif from the coasts
    else:
        g = grp

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

    # set boundary

    # open (water) boundaries

    try:
        water = b.boundary.geoms[0] - (b.boundary.geoms[0] - grl)
    except:
        water = b.boundary - (b.boundary - grl)

    try:
        cwater = shapely.ops.linemerge(water)
    except:
        cwater = water

    df_water = gp.GeoDataFrame(geometry=[cwater]).explode(index_parts=True).droplevel(0).reset_index(drop=True)
    df_water["tag"] = "open"

    # land boundaries!!
    try:
        land = b.boundary.geoms[0] - grl
    except:
        land = b.boundary - grl

    try:
        cland = shapely.ops.linemerge(land)
    except:
        cland = land

    df_land = gp.GeoDataFrame(geometry=[cland]).explode(index_parts=True).droplevel(0).reset_index(drop=True)
    df_land["tag"] = "land"

    ### interpolate if needed

    # even = kwargs.get("even", False)
    # ds = kwargs.get("ds", 0.001)

    #    if even:

    #        ibs = df_land.index.levels[0].shape[0]

    #        if ibs < 0:
    #            ## resample boundaries
    # ndfs = {}
    # for ic in tqdm(range(ibs)):
    #     contour = df_land.index.levels[0][ic]
    #     curve = df_land.loc[contour, ["lon", "lat"]]
    #     di = spline(curve, ds=ds)
    #     di["z"] = df_land.loc[contour].z.values[0]
    #     di["tag"] = df_land.loc[contour].tag.values[0].astype(int)
    #     ndfs.update({contour: di})
    #
    # df_land = pd.concat(ndfs, axis=0)
    # df_land["z"] = df_land.z.values.astype(int)
    # df_land["tag"] = df_land.tag.values.astype(int)

    ndf = (
        gp.GeoDataFrame(geometry=[shapely.ops.LineString(x) for x in b.interiors])
        .explode(index_parts=True)
        .droplevel(0)
        .reset_index(drop=True)
    )
    ndf["tag"] = "island"

    # ### interpolate if needed
    # if even:
    #
    #     ibs = ndf.index.levels[0].shape[0]
    #
    #     ## resample boundaries
    #     ndfs = {}
    #     for ic in tqdm(range(ibs)):
    #         contour = ndf.index.levels[0][ic]
    #         curve = ndf.loc[contour, ["lon", "lat"]]
    #         curve = pd.concat([curve, curve.loc[0:0]]).reset_index(
    #             drop=True
    #         )  # add the first point to do a better interpolation
    #         di = spline(curve, ds=ds)
    #         di["z"] = ndf.loc[contour].z.values[0]
    #         di["tag"] = ndf.loc[contour].tag.values[0].astype(int)
    #         ndfs.update({contour: di.drop_duplicates(["lon", "lat"])})  # remove duplicated points (see above)
    #
    #     ndf = pd.concat(ndfs, axis=0)
    #
    #     ndf["z"] = ndf.z.values.astype(int)
    #     ndf["tag"] = ndf.tag.values.astype(int)

    df = pd.concat([df_water, df_land, ndf]).reset_index(drop=True)

    df = df.loc[~df.is_empty]  # clean up

    return df


def global_tag(geo, cbuffer, blevels, R=1):

    # Manage coastlines
    logger.info("preparing coastlines")

    # buffer the coastlines on demand
    if cbuffer:
        logger.info("Buffering...")
        geo = buffer_(geo, cbuffer)
        logger.info("...done")

    # ANTARTICA
    anta_mask = geo.bounds.miny < geo.bounds.miny.min() + 0.1  # indentify antartica
    anta = geo.loc[anta_mask]
    indx = anta.index  # keep index

    if anta.geom_type.iloc[0] == "Polygon":  # convert boundary values to pandas
        anta = pd.DataFrame(anta.boundary.values[0].coords[:], columns=["lon", "lat"])
    else:
        anta = pd.DataFrame(anta.geometry.iloc[0].coords, columns=["lon", "lat"])

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
    m1 = geo[geo.bounds.minx == geo.bounds.minx.min()].index
    m2 = geo[geo.bounds.maxx == geo.bounds.maxx.max()].index
    mm = np.concatenate((m1, m2))  # join them
    mm = [j for j in mm if j != indx]  # subtract antartica

    # convert to u,v (stereographic coordinates)
    for idx, poly in geo.iterrows():
        if idx == indx:
            # resolve potential singularity
            gk = geo.loc[indx].geometry.iloc[0]
            gx = np.array([x for (x, y) in gk.coords[:]])
            gy = np.array([y for (x, y) in gk.coords[:]])
            jj = np.argwhere(gy > -86)
            gx_, gy_ = to_stereo(gx[jj], gy[jj], R=1)  # project valid values only
            geo.loc[idx, "geometry"] = shapely.geometry.LineString(list(zip(gx_, gy_)))
        else:
            geo.loc[idx, "geometry"] = shapely.ops.transform(
                lambda x, y, z=None: to_stereo(x, y, R=R),
                poly.geometry,
            )

    w = geo.drop(indx)  # get all polygons

    # join the split polygons
    ww = w.loc[mm]  # split entities
    qq = shapely.ops.polygonize_full(ww.geometry.values)  # polygonize in case of LineStrings
    if len(qq[0]) > 0:
        ww = gp.GeoDataFrame(geometry=list(qq[0]))  # convert to gp
    elif len(qq[2]) > 0:
        ww = gp.GeoDataFrame(geometry=list(qq[2]))  # convert to gp

    cs = []  # adjust values around zero (in projection - international line in Platee Carree)
    for idx, line in ww.itertuples():
        try:
            x_ = [x for (x, y) in line.exterior.coords[:]]  # polygons
            y_ = [y if abs(y) > 1.0e-4 else 0 for (x, y) in line.exterior.coords[:]]
        except:
            x_ = [x for (x, y) in line.coords[:]]  # LineStrings
            y_ = [y if abs(y) > 1.0e-4 else 0 for (x, y) in line.coords[:]]

        cs.append(shapely.geometry.Polygon(list(zip(x_, y_))))

    ww = gp.GeoDataFrame(geometry=cs)

    gw = gp.GeoDataFrame(
        geometry=list(ww.buffer(0).unary_union)
    )  # merge the polygons that are split (around -180/180)

    gw = gp.GeoDataFrame(geometry=gw.boundary.values)

    w = w.drop(mm)

    wp = w.loc[w.geom_type == "Polygon"]
    w2 = w.loc[w.geom_type != "Polygon"]
    w1 = gp.GeoDataFrame(geometry=wp.exterior.values)  # get boundaries if Polygons

    # Check antartica LineString
    if not geo.iloc[indx].geometry.values[0].is_ring:
        ca = gp.GeoDataFrame(
            geometry=[shapely.geometry.LinearRing(geo.loc[indx].geometry.values[0])],
            index=indx,
        )
        ca["geometry"] = shapely.geometry.LineString(ca.geometry.values[0])
    else:
        ca = geo.loc[indx]

    # PUT ALL TOGETHER
    geo = pd.concat([w1, w2, gw, ca], ignore_index=True).reset_index(drop=True)

    logger.info("storing boundaries")

    geo["tag"] = "island"

    geo["length"] = geo["geometry"][:].length  # sort with length
    geo = geo.sort_values(by="length", ascending=0)

    geo = gp.GeoDataFrame(geo)  # clean up

    # idx = 0
    # dic = {}
    # for i, line in tqdm(geo.iterrows(), total=geo.shape[0]):
    #     lon = []
    #     lat = []
    #
    #     for x, y in line.geometry.coords[:]:
    #         lon.append(x)
    #         lat.append(y)
    #     dic.update({"line{}".format(idx): {"lon": lon, "lat": lat, "tag": line.tag}})
    #     idx += 1
    #
    # # Serialize data into file:
    # json.dump( dic, open( "dich.json", 'w' ) )
    # sys.exit()
    #
    # dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}
    #
    # df = pd.concat(dict_of_df, axis=0)
    #
    # df["z"] = 0
    # df = df.drop_duplicates()  # drop the repeat value on closed boundaries

    return geo.loc[~geo.is_empty]
