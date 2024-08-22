import pyposeidon.dem as pdem
import pytest
import geopandas as gp
import cartopy.feature as cf


from . import DATA_DIR

COAST_FILE = (DATA_DIR / "ocean.parquet").as_posix()

DEM_SOURCES = pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
    ],
)


@pytest.fixture(scope="session")
def coasts():
    #    coast = gp.read_file(COAST_FILE)
    cr = "h"
    coast = cf.NaturalEarthFeature(
        category="physical", name="land", scale="{}m".format({"l": 110, "i": 50, "h": 10}[cr])
    )
    coast = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])
    return coast


@pytest.mark.slow
@DEM_SOURCES
def test_dem_adjust(coasts, dem_source):
    # Just elevation
    df1 = pdem.Dem(dem_source=dem_source)  # get dem
    check1 = df1.adjust(coasts, tiles=False)
    df2 = pdem.Dem(dem_source=dem_source)
    check2 = df2.adjust(coasts, tiles=True)

    assert check1 == check2 & df1.Dataset.equals(df2.Dataset)
