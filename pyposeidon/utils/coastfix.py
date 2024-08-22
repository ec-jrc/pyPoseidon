import sys
import shapely
import geopandas as gp

# logging setup
import logging

logger = logging.getLogger(__name__)


def simplify(geo):
    logger.info("simplify coastlines dataset")
    geo = geo.loc[:, ["geometry"]].explode(index_parts=True).droplevel(0).reset_index(drop=True)

    # Just in case merge
    if (geo.geom_type == "LineString").all():
        geo_ = gp.GeoDataFrame(geometry=[shapely.ops.linemerge(list(geo.geometry.values))])
        geo = geo_.explode(index_parts=True).droplevel(0).reset_index(drop=True)

        if geo.is_valid.all():
            return geo
        else:
            logger.error("Coastlines not valid... exiting")
            sys.exit(1)

    if (geo.geom_type == "Polygon").all():
        try:
            geo_ = list(geo.buffer(0).union_all().geoms)
        except TypeError:
            geo_ = [geo.buffer(0).union_all()]

        geo = gp.GeoDataFrame(geometry=geo_)
        geo = gp.GeoDataFrame(geometry=geo.buffer(0))

        if geo.is_valid.all() and (geo.boundary.geom_type == "LineString").all():
            return geo

        if (geo.boundary.geom_type == "MultiLineString").any():
            dg = geo.loc[geo.boundary.geom_type == "MultiLineString"]
            for idx, geom in dg.iterrows():
                pl = shapely.polygonize_full(dg.loc[[idx]].boundary.explode(index_parts=False).values)
                df = gp.GeoDataFrame(geometry=[pl[0]]).explode(index_parts=False)
                dff = gp.GeoDataFrame(geometry=[df.union_all()])
                geo.loc[idx, "geometry"] = dff.geometry.buffer(0).values[0]

            if geo.is_valid.all() and (geo.boundary.geom_type == "LineString").all():
                return geo
            else:
                logger.error("Coastlines not valid... exiting")
                sys.exit(1)
