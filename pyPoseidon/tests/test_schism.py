import pytest
import pyPoseidon
import os

PWD = os.getcwd()


#define in a dictionary the properties of the model..
case1={'solver':'schism',
     'grid_file': PWD + '/data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2017-10-1 0:0:0',
     'time_frame':'12H',
     'meteo_source' : [PWD + '/data/erai.grib'], #meteo file
     'dem_source' : PWD + '/data/dem.nc',
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

case2={'solver':'schism',
     'grid_file': PWD + '/data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2017-10-1 0:0:0',
     'end_date':'2017-10-2 0:0:0', #optional instead of time_frame
     'dem_source' : PWD + '/data/dem.nc',
     'meteo_source' : [PWD + '/data/erai.grib'], #meteo file
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

    
case3={'solver':'schism',
     'grid_file': PWD + '/data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2011-1-1 0:0:0',
     'time_frame':'12H',
     'meteo_source' : [PWD + '/data/era5.grib'], #meteo file
     'dem_source' : PWD + '/data/dem.nc',
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

def schism(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model
        
    b = pyPoseidon.model(**dic)
    
    try:
        b.execute()
        a = pyPoseidon.read_model(rpath+'test_model.json') # read model
        a.execute()
        return True
    except:
        return False

@pytest.mark.parametrize('case', [case1, case2, case3])
def test_answer(tmpdir, case):
    assert schism(tmpdir,case) == True