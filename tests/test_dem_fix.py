import pyposeidon.dem as pdem
import pyposeidon.mesh as pmesh
import pytest
import numpy as np
import geopandas as gp
import cartopy.feature as cf

from . import DATA_DIR

DEM_SOURCES = pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
        pytest.param(DATA_DIR / "dem.tif", id="local geotiff"),
    ],
)

WINDOWS = pytest.mark.parametrize(
    "window",
    [
        {"lon_min": -30, "lon_max": -10.0, "lat_min": 60.0, "lat_max": 70.0},
        {"lon_min": 175.0, "lon_max": 185.0, "lat_min": -21.5, "lat_max": -14.5},
        {"lon_min": -185.0, "lon_max": -175.0, "lat_min": -21.5, "lat_max": -14.5},
        {"lon_min": -3, "lon_max": 1.0, "lat_min": 42.0, "lat_max": 45.0},
    ],
)


@pytest.fixture(scope="session")
def natural_earth():
    coast = cf.NaturalEarthFeature(category="physical", name="land", scale="110m")
    gdf = gp.GeoDataFrame(geometry=[x.buffer(0.001) for x in coast.geometries()])
    return gdf


@DEM_SOURCES
@WINDOWS
def test_dem_adjust(natural_earth, dem_source, window):
    # Just elevation
    df = pdem.dem(**window, dem_source=dem_source)  # get dem
    df.adjust(natural_earth)
    assert np.isnan(df.Dataset.adjusted.values).sum() == 0


# Schism mesh
@pytest.mark.schism
@DEM_SOURCES
@WINDOWS
def test_schism_mesh(tmpdir, natural_earth, dem_source, window):
    mesh = pmesh.set(
        type="tri2d",
        geometry=window,
        coastlines=natural_earth,
        mesh_generator="jigsaw",
        rpath=str(tmpdir) + "/",
    )
    xg = mesh.Dataset.SCHISM_hgrid_node_x.values
    yg = mesh.Dataset.SCHISM_hgrid_node_y.values
    dem = pdem.dem(**window, dem_source=dem_source, adjust_dem=False)  # get dem
    dem.Dataset = pdem.dem_on_mesh(dem.Dataset, grid_x=xg, grid_y=yg)  # get dem on mesh
    dem.adjust(natural_earth)
    assert np.isnan(dem.Dataset.fval.values).sum() == 0


# D3D mesh
@pytest.mark.delft
@DEM_SOURCES
@WINDOWS
def test_d3d_mesh(tmpdir, natural_earth, dem_source, window):
    mesh = pmesh.set(type="r2d", geometry=window, resolution=0.1, rpath=str(tmpdir) + "/")
    gr = mesh.Dataset
    xp, yp = gr.lons.values, gr.lats.values
    dem = pdem.dem(**window, dem_source=dem_source, adjust_dem=False)  # get dem
    dem.Dataset = pdem.dem_on_mesh(dem.Dataset, grid_x=xp, grid_y=yp)  # get dem on mesh
    dem.adjust(natural_earth)
    assert np.isnan(dem.Dataset.fval.values).sum() == 0
