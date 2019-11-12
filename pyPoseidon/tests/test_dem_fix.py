import pyPoseidon.dem as pdem
import pytest
import numpy as np
import os

PWD = os.getcwd()


#define the lat/lon window and time frame of interest
window1 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.,
    'dem_source' : PWD + '/data/dem.nc'
}

window2={
     'lon_min':175., # lat/lon window
     'lon_max':185.,
     'lat_min':14.5,
     'lat_max':16.5,
     'dem_source' : PWD + '/data/dem.nc'
}

@pytest.mark.parametrize('kwargs', [ window1 , window2 ])
def test_answer(tmpdir, kwargs):
    ## lat,lon grid
    resolution=.1
    lon=np.arange(kwargs['lon_min'],kwargs['lon_max'],resolution)
    lat=np.arange(kwargs['lat_min'],kwargs['lat_max'],resolution)
    xp, yp = np.meshgrid(lon,lat)
    
    kwargs.update({'grid_x':xp, 'grid_y':yp})
    
    #get dem 
    df = pdem.dem(**kwargs)
    
    df.adjust(shpfile='./data/coast.shp',nc=50)
    
    assert np.isnan(df.Dataset.ival.values).sum() == 0
