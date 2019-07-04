import pytest
import pyPoseidon.model as pm

#define in a dictionary the properties of the model..
case1={'minlon' : -30,
     'maxlon' : -10.,
     'minlat' : 60.,
     'maxlat' : 70.,
     'start_date':'2018-10-1',
     'time_frame':'12H',
     'solver':'d3d',
     'resolution':0.2, #grid resoltuion 
     'map_step':20, # step for output of map field in d3d 
     'restart_step':60, # when to output restart file
     'ncores': 4 , #number of cores
     'meteo_files' : ['./data/oper/uvp_2018100100.grib'],
     'epath':'/Users/brey/DELFT3D/SVN/7545/bin/lnx64/', #folder for solver executables
     'conda_env':'pyPoseidon' # optional conda env for running the solver
#     'update':['all'] # optional to select update quantities
    }


def d3d(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model
    b = pm(**dic)

    try:
        b.execute()
        a = pm.read_model(rpath+'d3d_model.json') # read model
        a.execute()
        return True
    except:
        return False

@pytest.mark.parametrize('case', [case1])
def test_answer(tmpdir, case):
    assert d3d(tmpdir,case) == True