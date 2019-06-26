import pyPoseidon.dem as pdem
import numpy as np
import pytest

#define the lat/lon window and time frame of interest
window1 = {
    'minlon' : -30,
    'maxlon' : -10.,
    'minlat' : 60.,
    'maxlat' : 70.
}

window2={'minlon':175., # lat/lon window
     'maxlon':184.,
     'minlat':14.5,
     'maxlat':16.5}

window3={'minlon':175.-360., # lat/lon window
     'maxlon':184.-360.,
     'minlat':14.5,
     'maxlat':16.5}

@pytest.mark.parametrize('kwargs', [ window1, window2, window3 ])
def test_answer(tmpdir, kwargs):
    
    df = pdem.dem(**kwargs)
    
    check = np.isnan(df.altimetry.values).sum() == 0
    assert check == True