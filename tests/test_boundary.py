import pyposeidon.boundary as pb
import pytest
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

boundary = DATA_DIR / "bl"

land = cf.NaturalEarthFeature(category="physical", name="land", scale="50m")
nland = gp.GeoDataFrame(geometry=[x for x in land.geometries()])

coast = cf.NaturalEarthFeature(category="physical", name="coastline", scale="50m")
ncoast = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])

#gh = cf.GSHHSFeature(scale="intermediate", levels=[1])
#iGSHHS = gp.GeoDataFrame(geometry=[x for x in gh.geometries()])

INPUTS = pytest.mark.parametrize("input", [nland, ncoast]) #, iGSHHS])
WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 175.0, "lon_max": 184.0, "lat_min": 14.5, "lat_max": 16.5},
    ],
)

@INPUTS
@WINDOWS
def test_window(tmpdir, window, input):

    df = pb.get_boundaries(geometry=window, coastlines=input, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)

@INPUTS
def test_global(tmpdir, input):

    df = pb.get_boundaries(geometry='global', coastlines=input, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


@INPUTS
@WINDOWS
def test_buffer(tmpdir, window, input):

    df = pb.get_boundaries(geometry=window, coastlines=input, cbuffer = 0.01, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


@INPUTS
@WINDOWS
def test_isodem(tmpdir, window, input):

    df = pb.get_boundaries(geometry=window, blevels = [-100], rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


def test_custom(tmpdir):
    
    df = pb.get_boundaries(geometry=boundary, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)
    