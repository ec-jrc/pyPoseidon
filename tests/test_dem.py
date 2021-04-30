import pyposeidon.dem as pdem
import numpy as np
import pytest

from . import DATA_DIR


# define the lat/lon window and time frame of interest


@pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
        pytest.param(DATA_DIR / "dem.tif", id="local geotiff"),
        pytest.param("https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm30plus", id="url"),
    ],
)
@pytest.mark.parametrize(
    "kwargs",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 175.0, "lon_max": 184.0, "lat_min": 14.5, "lat_max": 16.5},
        {"lon_min": 175.0 - 360.0, "lon_max": 184.0 - 360.0, "lat_min": 14.5, "lat_max": 16.5},
    ],
)
def test_answer(dem_source, kwargs):
    df = pdem.dem(dem_source=dem_source, **kwargs)
    assert np.isnan(df.Dataset.elevation.values).sum() == 0
