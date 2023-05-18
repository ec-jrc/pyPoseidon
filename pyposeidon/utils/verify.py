#!/usr/bin/env python
# coding: utf-8

import os

import numpy as np
import pandas as pd
import geopandas as gp
import xarray as xr
from tqdm import tqdm
import shapely

import logging

logger = logging.getLogger(__name__)


def verify(g, shp, thorough=False):
    # ---------------------------------------------------------------------
    logger.info(" Verify mesh against coastline\n")
    # ---------------------------------------------------------------------

    lon_min = g.Dataset.SCHISM_hgrid_node_x.values.min()
    lon_max = g.Dataset.SCHISM_hgrid_node_x.values.max()
    lat_min = g.Dataset.SCHISM_hgrid_node_y.values.min()
    lat_max = g.Dataset.SCHISM_hgrid_node_y.values.max()

    c = shp.cx[lon_min:lon_max, lat_min:lat_max]

    # ## Test polygons

    d = g.Dataset

    x = d.SCHISM_hgrid_node_x.values
    y = d.SCHISM_hgrid_node_y.values
    tri = d.SCHISM_hgrid_face_nodes.values

    nodes = pd.DataFrame({"lon": x, "lat": y})

    elems = pd.DataFrame(tri, columns=["a", "b", "c"])

    bnodes = g.Dataset[["node", "id", "type"]].to_dataframe()

    # ### Find the invalid nodes (that cross the coasts)
    cos = c.geometry
    cos_ = shapely.set_operations.union_all(cos)

    gps = shapely.points(list(nodes.values))

    gtree = shapely.STRtree(gps)

    invs = gtree.query(cos_, predicate="contains").tolist()

    # ---------------------------------------------------------------------
    logger.info("Number of nodes within/touching the coastlines {}\n".format(len(invs)))
    # ---------------------------------------------------------------------

    nps = len(invs)

    nels = 1

    if thorough:
        # ### Find invalid elements (that cross land)

        # cells to polygons
        ap = nodes.loc[elems.a]
        bp = nodes.loc[elems.b]
        cp = nodes.loc[elems.c]

        elems["ap"] = ap.values.tolist()
        elems["bp"] = bp.values.tolist()
        elems["cp"] = cp.values.tolist()

        n = 2
        al = elems.ap + elems.bp + elems.cp + elems.ap
        coords = [[l[i : i + n] for i in range(0, len(l), n)] for l in al]
        elems["coordinates"] = coords

        jig = shapely.polygons(coords)

        jtree = shapely.STRtree(jig)

        jig_ = shapely.set_operations.union_all(jig)

        cross = shapely.set_operations.intersection(jig_, cos_)

        # #### convert to dataframe

        fd = gp.GeoDataFrame({"geometry": cross}, index=[0])

        ipols = fd.explode(index_parts=True).loc[0]

        # #### Reject small injuctions
        mask = ipols.area.values == 0.0

        ipols = ipols[~mask].reset_index(drop=True)
        ipols = gp.GeoDataFrame(ipols)

        # ---------------------------------------------------------------------
        logger.info("Number of elements intersecting the coastlines {}\n".format(ipols.shape[0]))
        # ---------------------------------------------------------------------
        logger.info("Maximum intersecting area {}\n".format(ipols.area.max()))
        # ---------------------------------------------------------------------

        nels = ipols.shape[0]

        if ipols.area.max() < 1e-3:
            # ---------------------------------------------------------------------
            logger.info("Mesh is verified against the coastline")
            # ---------------------------------------------------------------------
            return True

    if nps == 0 and nels == 0:
        # ---------------------------------------------------------------------
        logger.info("Mesh is verified against the coastline")
        # ---------------------------------------------------------------------
        return True
    elif nps == 0:
        # ---------------------------------------------------------------------
        logger.info("Mesh is node verified against the coastline")
        # ---------------------------------------------------------------------
        return True
    else:
        # ---------------------------------------------------------------------
        logger.warning("Mesh can't be verified against the coastline")
        # ---------------------------------------------------------------------
        return False
