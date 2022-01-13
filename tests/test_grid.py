import pyposeidon.mesh as pmesh
import numpy as np
import pytest
import os
import geopandas as gp
import cartopy.feature as cf


cr = "i"
coast = cf.NaturalEarthFeature(
    category="physical",
    name="land",
    scale="{}m".format({"l": 110, "i": 50, "h": 10}[cr]),
)


natural_earth = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])

coast = cf.GSHHSFeature(scale="auto", levels=[1])

GSHHS = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])


# define the lat/lon window and time frame of interest
window0 = {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0}

window1 = {
    "lon_min": 175.0,
    "lon_max": 184.0,
    "lat_min": -21.5,
    "lat_max": -14.5,
}  # lat/lon window

window2 = {"lon_min": -20.0, "lon_max": -10.0, "lat_min": 63.0, "lat_max": 67.0}

window3 = {
    "lon_min": 175.0 - 360.0,
    "lon_max": 184.0 - 360.0,
    "lat_min": -21.5,
    "lat_max": -14.5,
}  # lat/lon window

window4 = {"lon_min": -25.0, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 68.0}


@pytest.mark.parametrize("window", [window0, window1, window2, window3, window4])
@pytest.mark.parametrize("coast", [natural_earth, GSHHS])
@pytest.mark.parametrize("ggor", ["jigsaw", "gmsh"])
def test_answer(tmpdir, window, coast, ggor):

    mesh = pmesh.set(
        type="tri2d",
        geometry=window,
        coastlines=coast,
        rpath=str(tmpdir) + "/",
        mesh_generator=ggor,
    )

    check = np.isnan(mesh.Dataset.depth.values).sum() == 0
    assert check == True
