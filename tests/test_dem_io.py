import pyposeidon.dem as pdem
import pyposeidon.model as pmodel
import pyposeidon.mesh as pmesh


import pytest
import xarray as xr
import os
import shutil
import numpy as np


from . import DATA_DIR

DEM_SOURCES = pytest.mark.parametrize(
    "dem_source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="local netcdf"),
        pytest.param(DATA_DIR / "dem.tif", id="local geotiff"),
    ],
)


# define the lat/lon window and time frame of interest
window1 = {
    "lon_min": -30,
    "lon_max": -10.0,
    "lat_min": 60.0,
    "lat_max": 70.0,
}


def schism(tmpdir, dem_source, kwargs):

    mesh = pmesh.set(type="tri2d", mesh_file=(DATA_DIR / "hgrid.gr3").as_posix())

    # update kwargs
    xp = mesh.Dataset.SCHISM_hgrid_node_x.values
    yp = mesh.Dataset.SCHISM_hgrid_node_y.values

    kwargs.update({"grid_x": xp, "grid_y": yp})

    # get dem
    df = pdem.dem(dem_source=dem_source, **kwargs)

    # get dem on mesh
    df.Dataset = pdem.dem_on_mesh(df.Dataset, **kwargs)

    mesh.Dataset["depth"].loc[:] = -df.Dataset.ival.values

    filename_ = str(tmpdir.join("hgrid_.gr3"))
    # output to mesh file
    mesh.to_file(filename_)

    # read again new mesh
    mesh_ = pmesh.set(type="tri2d", mesh_file=filename_)

    # compare
    return mesh.Dataset.equals(mesh_.Dataset)


def d3d(tmpdir, dem_source, kwargs):

    ## lat,lon grid
    resolution = 0.1
    lon = np.arange(kwargs["lon_min"], kwargs["lon_max"], resolution)
    lat = np.arange(kwargs["lat_min"], kwargs["lat_max"], resolution)
    xp, yp = np.meshgrid(lon, lat)
    # update kwargs
    kwargs.update({"grid_x": xp, "grid_y": yp})

    # get dem
    df = pdem.dem(dem_source=dem_source, **kwargs)
    # get dem on mesh
    df.Dataset = pdem.dem_on_mesh(df.Dataset, **kwargs)

    rpath = str(tmpdir) + "/"
    # output
    pdem.to_output(df.Dataset, solver="d3d", rpath=rpath)

    # read again dem
    m = pmodel.set(solver="d3d")
    rd = m.from_dep(rpath + "d3d.dep")

    # compare
    c1 = -rd.where(rd != -999)
    c2 = df.Dataset.ival.where(df.Dataset.ival < 0)

    return c1.fillna(0).equals(c2.fillna(0))


@DEM_SOURCES
@pytest.mark.parametrize("kwargs", [window1])
@pytest.mark.parametrize("solver", [schism, d3d])
def test_answer(tmpdir, dem_source, kwargs, solver):
    assert solver(tmpdir, dem_source, kwargs) == True
