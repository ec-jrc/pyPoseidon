import pyPoseidon.meteo as pmeteo
import pytest
import xarray as xr
import os
import shutil
import numpy as np


def func(tmpdir,name):

    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', name)
    #read meteo file
    df = pmeteo.cfgrib(filenames=filename)
    df = df.sortby('latitude', ascending=True)

    rpath = str(tmpdir)+'/'
    #output to uvp files
    pmeteo.to_output(df,solver='schism',rpath=rpath)

    #read again meteo
    path = rpath + '/sflux/'
    dr = xr.open_dataset(path + '/sflux_air_1.001.nc')

    #cleanup
    try:
        shutil.rmtree(path)
    except OSError as e:
        print ("Error: %s - %s." % (e.filename, e.strerror))

    #compare
    msl = np.array_equal(df.msl.values,dr.prmsl.values)
    u10 = np.array_equal(df.u10.values,dr.uwind.values)
    v10 = np.array_equal(df.v10.values,dr.vwind.values)

    return all([msl,u10,v10])
    
@pytest.mark.parametrize('filename', ['erai.grib', 'era5.grib'])
def test_answer(tmpdir, filename):
    assert func(tmpdir,filename) == True
