import pyposeidon.boundary as pb
import pytest
import geopandas as gp

from . import DATA_DIR

noaa = DATA_DIR / "bl.zip"

COAST_FILE = (DATA_DIR / "ocean.parquet").as_posix()

land = gp.read_file(COAST_FILE)
coast = gp.GeoDataFrame(geometry=land.boundary)

INPUTS = pytest.mark.parametrize("input", [land, coast])
CUSTOM = pytest.mark.parametrize("boundary", [noaa])
WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 175.0, "lon_max": 184.0, "lat_min": -50, "lat_max": -30},
    ],
)
DEM_SOURCES = pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
        pytest.param(DATA_DIR / "dem.tif", id="local geotiff"),
    ],
)


@INPUTS
@WINDOWS
def test_window(tmpdir, window, input):
    df = pb.Boundary(geometry=window, coastlines=input, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


@INPUTS
def test_global(tmpdir, input):
    df = pb.Boundary(geometry="global", coastlines=input, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


@INPUTS
@WINDOWS
def test_buffer(tmpdir, window, input):
    df = pb.Boundary(geometry=window, coastlines=input, cbuffer=0.01, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)


# @INPUTS
# @WINDOWS
# @DEM_SOURCES
# def test_isodem(tmpdir, window, input, dem_source):
#    df = pb.Boundary(geometry=window, dem_source=dem_source, blevels=[-100], rpath=str(tmpdir) + "/")
#    assert isinstance(df.contours, gp.GeoDataFrame)


# @pytest.mark.slow
# def test_isodem_with_url(tmpdir):
#    input = nland
#    window = {"lon_min": 176.5, "lon_max": 177.0, "lat_min": 16.0, "lat_max": 16.5}
#    dem_source = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm30plus"
#    df = pb.Boundary(geometry=window, dem_source=dem_source, blevels=[-100], rpath=str(tmpdir) + "/")
#    assert isinstance(df.contours, gp.GeoDataFrame)


@CUSTOM
def test_custom(tmpdir, boundary):
    df = pb.Boundary(geometry=boundary, rpath=str(tmpdir) + "/")

    assert isinstance(df.contours, gp.GeoDataFrame)
