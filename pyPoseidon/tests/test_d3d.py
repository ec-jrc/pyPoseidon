import pytest
import pyPoseidon
import os

PWD = os.getcwd()

#define in a dictionary the properties of the model..
case1={'lon_min' : -30,
     'lon_max' : -10.,
     'lat_min' : 60.,
     'lat_max' : 70.,
     'start_date':'2018-10-1',
     'time_frame':'12H',
     'solver':'d3d',
     'resolution':0.2, #grid resoltuion 
     'map_step':20, # step for output of map field in d3d 
     'restart_step':60, # when to output restart file
     'ncores': 4 , #number of cores
     'meteo_source' : [PWD + '/data/uvp_2018100100.grib'],
     'dem_source' : PWD + '/data/dem.nc',
     'engine':'cfgrib',
#     'update':['all'] # optional to select update quantities
    }


def d3d(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model
    b = pyPoseidon.model(**dic)

    try:
        b.execute()
        a = pyPoseidon.read_model(rpath+'d3d_model.json') # read model
        a.execute()
        return True
    except:
        return False

@pytest.mark.parametrize('case', [case1])
def test_answer(tmpdir, case):
    assert d3d(tmpdir,case) == True