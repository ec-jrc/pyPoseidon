import pyPoseidon.meteo as pmeteo
import pytest
import xarray as xr
import os
import datetime
import numpy as np


kwargs = {
    'minlon' : -30,
    'maxlon' : -10.,
    'minlat' : 60.,
    'maxlat' : 70., 
}

def schism(tmpdir,df, kwargs):

    #get meteo 
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
    
def d3d(tmpdir,df, kwargs):

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
    
    
    
@pytest.mark.parametrize('url', [ 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best'] ) #, 'https://cidportal.jrc.ec.europa.eu/services/erddap/griddap/NCEP_Global_Best' ])
def test_answer(tmpdir, url): 
    start_date = datetime.datetime.today() +  datetime.timedelta(days=-30)  
    end_date = start_date + datetime.timedelta(days=1)
      
    start_date = start_date.strftime(format='%Y-%m-%d')
    end_date = end_date.strftime(format='%Y-%m-%d')
    kwargs.update({'start_date':start_date, 'end_date':end_date})

    #get meteo     
    df = pmeteo.from_url(url=url, **kwargs)
    

    assert schism(tmpdir,df,kwargs) == True
    assert d3d(tmpdir,df,kwargs) == True
    
