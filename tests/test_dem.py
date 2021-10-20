import pyposeidon.dem as pdem
import numpy as np
import pytest

from . import DATA_DIR


DEM_SOURCES = pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
        pytest.param(DATA_DIR / "dem.tif", id="local geotiff"),
    ],
)

WINDOWS = pytest.mark.parametrize(
    "kwargs",
    [
        # fmt: off
        pytest.param({"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0}, id="local 1"),
        pytest.param({"lon_min": 175.0, "lon_max": 184.0, "lat_min": 14.5, "lat_max": 16.5}, id="local 2"),
        pytest.param({"lon_min": 175.0 - 360.0, "lon_max": 184.0 - 360.0, "lat_min": 14.5, "lat_max": 16.5}, id="local 3"),
        pytest.param({"lon_min": -180.0, "lon_max": 180.0, "lat_min": -90.0, "lat_max": 90.0}, id="global 2"),
        pytest.param({"lon_min": -360.0, "lon_max": 0, "lat_min": -90.0, "lat_max": 90.0}, id="global 1"),
        # fmt: on
    ],
)


@DEM_SOURCES
@WINDOWS
def test_dem(dem_source, kwargs):
    df = pdem.dem(dem_source=dem_source, **kwargs)
    assert np.isnan(df.Dataset.elevation.values).sum() == 0


def test_dem_source_is_url():
    dem_source = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm30plus"
    kwargs = {"lon_min": 176.5, "lon_max": 177.0, "lat_min": 16.0, "lat_max": 16.5}
    df = pdem.dem(dem_source=dem_source, **kwargs)
    assert np.isnan(df.Dataset.elevation.values).sum() == 0


def test_dem_source_is_not_specified():
    kwargs = {"lon_min": 176.5, "lon_max": 177.0, "lat_min": 16.0, "lat_max": 16.5}
    df = pdem.dem(**kwargs)
    assert np.isnan(df.Dataset.elevation.values).sum() == 0
