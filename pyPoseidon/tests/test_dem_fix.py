import pyPoseidon.dem as pdem
import pytest
import numpy as np


#define the lat/lon window and time frame of interest
window1 = {
    'minlon' : -30,
    'maxlon' : -10.,
    'minlat' : 60.,
    'maxlat' : 70.
}

window2={'minlon':175., # lat/lon window
     'maxlon':185.,
     'minlat':14.5,
     'maxlat':16.5}

@pytest.mark.parametrize('kwargs', [ window1 , window2 ])
def test_answer(tmpdir, kwargs):
    ## lat,lon grid
    resolution=.1
    lon=np.arange(kwargs['minlon'],kwargs['maxlon'],resolution)
    lat=np.arange(kwargs['minlat'],kwargs['maxlat'],resolution)
    xp, yp = np.meshgrid(lon,lat)
    
    kwargs.update({'grid_x':xp, 'grid_y':yp})
    
    #get dem 
    df = pdem.dem(**kwargs)
    
    df.adjust(shpfile='./data/coast.shp',nc=50)
    
    assert np.isnan(df.altimetry.ival.values).sum() == 0
