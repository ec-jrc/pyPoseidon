import pytest
import pyPoseidon
import os

from . import DATA_DIR

DEM_FILE = DATA_DIR / 'dem.nc'

case0 = {
    'lon_min' : -30,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 70.}

case1={'lon_min':175., # lat/lon window
     'lon_max':184.,
     'lat_min':-21.5,
     'lat_max':-14.5}
     
case2 = {
     'lon_min' : -20.,
     'lon_max' : -10.,
     'lat_min' : 63.,
     'lat_max' : 67.}

case3={'lon_min':175.-360., # lat/lon window
     'lon_max':184.-360.,
     'lat_min':-21.5,
     'lat_max':-14.5}
     
case4 = {
    'lon_min' : -25.,
    'lon_max' : -10.,
    'lat_min' : 60.,
    'lat_max' : 68.}


def schism(tmpdir,case):
    #initialize a model
    dic={'solver':'schism',
         'geometry': case,
         'manning':.12,
         'windrot':0.00001,
         'tag':'test',
         'start_date':'2017-10-1 0:0:0',
         'time_frame':'12H',
         'meteo_source' : [DATA_DIR / 'erai.grib'], #meteo file
         'engine':'cfgrib',
         'dem_source' : DEM_FILE,
         'ncores': 4 , #number of cores
         'update':['all'], #update only meteo, keep dem
         'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
        }
        
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model

    b = pyPoseidon.model(**dic)

    try:
        b.execute()
        return True
    except:
        return False

@pytest.mark.solvers
@pytest.mark.parametrize('case', [case0, case1, case2, case3])
def test_answer(tmpdir, case):
    assert schism(tmpdir,case) == True
