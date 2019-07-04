import pytest
import pyPoseidon.model as pm

#define in a dictionary the properties of the model..
case1={'solver':'schism',
     'grid_file': './data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2017-1-1 0:0:0',
     'time_frame':'12H',
     'ncores': 4 , #number of cores
     'conda_env':'pyPoseidon',
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

case2={'solver':'schism',
     'grid_file': './data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2017-1-1 0:0:0',
     'end_date':'2017-1-2 0:0:0', #optional instead of time_frame
     'ncores': 4 , #number of cores
     'conda_env':'pyPoseidon',
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

    
case3={'solver':'schism',
     'grid_file': './data/hgrid.gr3', 
     'manning':.12,
     'windrot':0.00001,
     'tag':'test',
     'start_date':'2011-1-1 0:0:0',
     'time_frame':'12H',
     'meteo_files' : ['./data/era5.grib'], #meteo file
     'ncores': 4 , #number of cores
     'conda_env':'pyPoseidon',
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':0.3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

def schism(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model
        
    b = pm(**dic)
    
    try:
        b.execute()
        a = pm.read_model(rpath+'test_model.json') # read model
        a.execute()
        return True
    except:
        return False

@pytest.mark.parametrize('case', [case1, case2, case3])
def test_answer(tmpdir, case):
    assert schism(tmpdir,case) == True