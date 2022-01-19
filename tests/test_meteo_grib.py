import pyposeidon.meteo as pmeteo
import pyposeidon.model as pmodel
import pytest
import xarray as xr
import os
import shutil
import numpy as np

from . import DATA_DIR


@pytest.mark.parametrize("input_name", ["erai.grib", "era5.grib", "uvp_2018100112.grib"])
def test_schism(tmpdir, input_name):
    filename = (DATA_DIR / input_name).as_posix()
    # read meteo file
    df = pmeteo.Meteo(meteo_source=filename)
    df.Dataset = df.Dataset.sortby("latitude", ascending=True)

    rpath = str(tmpdir) + "/"
    # output to uvp files
    df.to_output(solver_name="schism", rpath=rpath)

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


@pytest.mark.parametrize("input_name", ["erai.grib", "era5.grib", "uvp_2018100112.grib"])
def test_d3d(tmpdir, input_name):
    filename = (DATA_DIR / input_name).as_posix()
    # read meteo file
    df = pmeteo.Meteo(meteo_source=filename)

    rpath = str(tmpdir) + "/"
    # output to uvp files
    df.to_output(solver_name="d3d", rpath=rpath)

    # read again meteo
    m = pmodel.set(solver_name="d3d")

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
