import pyPoseidon.meteo as pmeteo
import pyPoseidon.model as pmodel
import pytest
import xarray as xr
import os
import shutil
import numpy as np


def schism(tmpdir,name):

    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', name)
    #read meteo file
    df = pmeteo.meteo(meteo_source=filename, engine='cfgrib')
    df.Dataset = df.Dataset.sortby('latitude', ascending=True)

    rpath = str(tmpdir)+'/'
    #output to uvp files
    df.to_output(solver='schism',rpath=rpath)

    #read again meteo
    path = rpath + '/sflux/'
    dr = xr.open_dataset(path + '/sflux_air_1.001.nc')

    #cleanup
#    try:
#        shutil.rmtree(path)
#    except OSError as e:
#        print ("Error: %s - %s." % (e.filename, e.strerror))

    #compare
    msl = np.array_equal(df.Dataset.msl.values,dr.prmsl.values)
    u10 = np.array_equal(df.Dataset.u10.values,dr.uwind.values)
    v10 = np.array_equal(df.Dataset.v10.values,dr.vwind.values)

    return all([msl,u10,v10])
    
def d3d(tmpdir,name):

    filename = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data', name)
    #read meteo file
    df = pmeteo.meteo(meteo_source=filename, engine='cfgrib')
    
    rpath = str(tmpdir)+'/'
    #output to uvp files
    df.to_output(solver='d3d',rpath=rpath)

    #read again meteo
    m = pmodel(solver='d3d')
    
    p = m.from_force(rpath + 'p.amp', 'msl')
    u = m.from_force(rpath + 'u.amu', 'u10')
    v = m.from_force(rpath + 'v.amv', 'v10')
    
    dr = xr.merge([p,u,v])
    dr = dr.sortby('latitude', ascending=True)
    
    #compare  
    df.Dataset = df.Dataset.sortby('latitude', ascending=True)
          
    msl = np.abs(df.Dataset.msl.values - dr.msl.values).max() < 1e-3
    u10 = np.abs(df.Dataset.u10.values - dr.u10.values).max() < 1e-3
    v10 = np.abs(df.Dataset.v10.values - dr.v10.values).max() < 1e-3

    return all([msl,u10,v10])
    
    
    
@pytest.mark.parametrize('filename', ['erai.grib', 'era5.grib', 'uvp_2018100112.grib' ])
def test_answer(tmpdir, filename):
    assert schism(tmpdir,filename) == True
    assert d3d(tmpdir,filename) == True
    
