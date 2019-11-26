import pyPoseidon.dem as pdem
import numpy as np
import pytest

from . import DATA_DIR

DEM_SOURCE = DATA_DIR / "dem.nc"


#define the lat/lon window and time frame of interest
window1 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.,
    'dem_source' : DEM_SOURCE
}

window2={'lon_min':175., # lat/lon window
     'lon_max':184.,
     'lat_min':14.5,
     'lat_max':16.5,
     'dem_source' : DEM_SOURCE,
}

window3={'lon_min':175.-360., # lat/lon window
     'lon_max':184.-360.,
     'lat_min':14.5,
     'lat_max':16.5,
     'dem_source' : DEM_SOURCE,
}

@pytest.mark.parametrize('kwargs', [ window1, window2, window3 ])
def test_answer(tmpdir, kwargs):
    
    df = pdem.dem(**kwargs)
    
    check = np.isnan(df.Dataset.elevation.values).sum() == 0
    assert check == True