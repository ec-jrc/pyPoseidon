import pyPoseidon.meteo as pmeteo
import pytest
import xarray as xr
import os
import shutil
import numpy as np


def schism(tmpdir,name):

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
#    try:
#        shutil.rmtree(path)
#    except OSError as e:
#        print ("Error: %s - %s." % (e.filename, e.strerror))

    #compare
    msl = np.array_equal(df.msl.values,dr.prmsl.values)
    u10 = np.array_equal(df.u10.values,dr.uwind.values)
    v10 = np.array_equal(df.v10.values,dr.vwind.values)

    return all([msl,u10,v10])
    
def d3d(tmpdir,name):

    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', name)
    #read meteo file
    df = pmeteo.cfgrib(filenames=filename)
    
    rpath = str(tmpdir)+'/'
    #output to uvp files
    pmeteo.to_output(df,solver='d3d',rpath=rpath)

    #read again meteo
    p = pmeteo.read_d3d_meteo(rpath + 'p.amp', 'msl')
    u = pmeteo.read_d3d_meteo(rpath + 'u.amu', 'u10')
    v = pmeteo.read_d3d_meteo(rpath + 'v.amv', 'v10')
    
    dr = xr.merge([p,u,v])
    dr = dr.sortby('latitude', ascending=True)
    
    #compare  
    df = df.sortby('latitude', ascending=True)
          
    msl = np.abs(df.msl.values - dr.msl.values).max() < 1e-3
    u10 = np.abs(df.u10.values - dr.u10.values).max() < 1e-3
    v10 = np.abs(df.v10.values - dr.v10.values).max() < 1e-3

    return all([msl,u10,v10])
    
    
    
@pytest.mark.parametrize('filename', ['erai.grib', 'era5.grib', 'oper/uvp_2018100112.grib' ])
def test_answer(tmpdir, filename):
    assert schism(tmpdir,filename) == True
    assert d3d(tmpdir,filename) == True
    
