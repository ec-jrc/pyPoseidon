import pyposeidon.meteo as pmeteo
import pyposeidon.model as pmodel
import pytest
import xarray as xr
import os
import shutil
import numpy as np

from . import DATA_DIR


def schism(tmpdir, name):
    filename = (DATA_DIR / name).as_posix()
    # read meteo file
    df = pmeteo.meteo(meteo_source=filename)
    df.Dataset = df.Dataset.sortby("latitude", ascending=True)

    rpath = str(tmpdir) + "/"
    # output to uvp files
    df.to_output(solver="schism", rpath=rpath)

    # read again meteo
    path = rpath + "/sflux/"
    dr = xr.open_dataset(path + "/sflux_air_1.0001.nc")

    # cleanup
    #    try:
    #        shutil.rmtree(path)
    #    except OSError as e:
    #        print ("Error: %s - %s." % (e.filename, e.strerror))

    # compare
    assert np.array_equal(df.Dataset.msl.values, dr.prmsl.values)
    assert np.array_equal(df.Dataset.u10.values, dr.uwind.values)
    assert np.array_equal(df.Dataset.v10.values, dr.vwind.values)


def d3d(tmpdir, name):
    filename = (DATA_DIR / name).as_posix()
    # read meteo file
    df = pmeteo.meteo(meteo_source=filename)

    rpath = str(tmpdir) + "/"
    # output to uvp files
    df.to_output(solver="d3d", rpath=rpath)

    # read again meteo
    m = pmodel.set(solver="d3d")

    p = m.from_force(rpath + "p.amp", "msl")
    u = m.from_force(rpath + "u.amu", "u10")
    v = m.from_force(rpath + "v.amv", "v10")

    dr = xr.merge([p, u, v])
    dr = dr.sortby("latitude", ascending=True)

    # compare
    df.Dataset = df.Dataset.sortby("latitude", ascending=True)

    assert np.abs(df.Dataset.msl.values - dr.msl.values).max() < 1e-3
    assert np.abs(df.Dataset.u10.values - dr.u10.values).max() < 1e-3
    assert np.abs(df.Dataset.v10.values - dr.v10.values).max() < 1e-3


@pytest.mark.parametrize("filename", ["erai.grib", "era5.grib", "uvp_2018100112.grib"])
@pytest.mark.parametrize("solver", [schism, d3d])
def test_meteo_grib(tmpdir, filename, solver):
    solver(tmpdir, filename)
