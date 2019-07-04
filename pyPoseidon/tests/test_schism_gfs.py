import pytest
import pyPoseidon.model as pm
import datetime
import numpy as np

#define in a dictionary the properties of the model..
case={'solver':'schism',
     'grid_file': '../pyPoseidon/tests/data/hgrid.gr3',
     'manning':.12,
     'windrot':0.00001,
     'tag':'iceland',
     'url':'https://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs{}/gfs_0p25_1hr_{}z'.format(start_date.strftime('%Y%m%d'),r[h]),
     'start_date': start_date.strftime(format='%Y-%m-%d %{}:0:0'.format(r[h])),
     'time_frame':'85H',
     'ncores': 4 , #number of cores
     'conda_env': 'pyPoseidon', # optional conda env for the solver
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':3, 'hotout':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'hotout_write':108 }
    }

def schism(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)+'/'
    dic.update({'rpath':rpath}) # use tmpdir for running the model
    
    #set the gfs run
    start_date = datetime.datetime.today() +  datetime.timedelta(days=-1)
    print(datetime.datetime.today())
    r = [0,6,12,18] # evaluate actual forecast hour available
    h = np.argmin([n for n in [start_date.hour-x for x in r]  if n>0])
    dic.update({'url':'https://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs{}/gfs_0p25_1hr_{}z'.format(start_date.strftime('%Y%m%d'),r[h])})

    b = pm(**dic)

    try:
        b.execute()
        return True
    except:
        return False

@pytest.mark.parametrize('case', [case])
def test_answer(tmpdir, case):
    assert schism(tmpdir,case) == True