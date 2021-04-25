# Create geojon features for animations
import pyposeidon
import pandas as pd
import geopandas as gp
import numpy as np
import branca
import shapely
import os
import json
import logging

logger = logging.getLogger("pyposeidon")


def create_geojson_features(gf, time):
    #    print('> Creating GeoJSON features...')
    features = []
    for idx, row in gf.iterrows():
        feature = {
            "id": idx,
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [[row["ap"].tolist(), row["bp"].tolist(), row["cp"].tolist()]],
            },
            "properties": {"times": [time], "style": {"fillColor": row["color"], "weight": 0, "fill-opacity": 1.0}},
        }
        features.append(feature)
    return features


def create_features(path="./", tag="schism"):

    logger.info("create geojson file for elevation")
    m = pyposeidon.read_model(path + tag + "_model.json")
    m.get_data()

    d = m.data.Dataset

    vmin = d.elev.min().compute()
    vmax = d.elev.max().compute()

    x = d.SCHISM_hgrid_node_x.values
    y = d.SCHISM_hgrid_node_y.values
    tri = d.SCHISM_hgrid_face_nodes.values

    nodes = pd.DataFrame({"lon": x, "lat": y})
    tria = pd.DataFrame(tri, columns=["a", "b", "c"])

    tria["ap"] = tria.apply(lambda x: nodes.loc[x.a, ["lon", "lat"]].values, axis=1)
    tria["bp"] = tria.apply(lambda x: nodes.loc[x.b, ["lon", "lat"]].values, axis=1)
    tria["cp"] = tria.apply(lambda x: nodes.loc[x.c, ["lon", "lat"]].values, axis=1)

    tria["geometry"] = tria.apply(lambda x: shapely.geometry.Polygon([x.ap, x.bp, x.cp]), axis=1)

    colormap = branca.colormap.LinearColormap(["green", "yellow", "red"], vmin=vmin.values, vmax=vmax.values)
    colormap.caption = "Elevation"

    # geopandas
    gf_ = gp.GeoDataFrame(tria, crs={"init": "epsg:4326"})

    gf_ = gf_.drop(["a", "b", "c", "ap", "bp", "cp"], axis=1)

    ## All frames
    fs = []
    for l in range(d.time.shape[0]):
        fr = d.elev.isel(time=l)
        a = fr.sel(nSCHISM_hgrid_node=tria.a.to_list()).values
        b = fr.sel(nSCHISM_hgrid_node=tria.b.to_list()).values
        c = fr.sel(nSCHISM_hgrid_node=tria.c.to_list()).values
        tria["value"] = np.mean([a, b, c], axis=0)
        tria["color"] = [colormap(x) for x in tria.value.to_list()]
        fs.append(create_geojson_features(tria, d.time[l].values.astype(str).split(".")[0]))

    tf = [j for i in fs for j in i]

    if not os.path.exists(path + "server"):
        os.makedirs(path + "server")

    json.dump(tf, open(path + "server/anim.json", "w"))

    gf_.to_file(path + "server/grid.geojson", driver="GeoJSON")

    logger.info("... saved")

    return
