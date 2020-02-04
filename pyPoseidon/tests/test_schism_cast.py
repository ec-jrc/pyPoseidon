import pytest
import pyPoseidon
from pyPoseidon.utils import cast, data
import json
import pandas as pd
import datetime
import os
import numpy as np

from . import DATA_DIR


GRID_FILE = (DATA_DIR / 'hgrid.gr3').as_posix()
DEM_FILE = (DATA_DIR / 'dem.nc').as_posix()
METEO_FILES_1 = [(DATA_DIR / name).as_posix() for name in ('uvp_2018100100.grib', 'uvp_2018100112.grib')]
METEO_FILES_2 = [(DATA_DIR / name).as_posix() for name in ('uvp_2018100100.grib', 'uvp_2018100112.grib', 'uvp_2018100200.grib', 'uvp_2018100212.grib')]


#define in a dictionary the properties of the model..
case={'solver':'schism',
     'grid_file': GRID_FILE,
     'manning':.12,
     'windrot':0.00001,
     'tag':'schism',
     'start_date':'2018-10-1 0:0:0',
     'time_frame':'12H',
     'dem_source' : DEM_FILE,
     'meteo_source' : METEO_FILES_1,
     'engine':'cfgrib',
     'combine_forecast': True, #combine meteo
     'combine_by':'nested',
     'xr_kwargs': {'concat_dim':'step'},
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':.5, 'nhot':1, 'ihot':0,'nspool':9, 'ihfskip':36, 'nhot_write':108 }
    }

#define in a dictionary the properties of the model..
check={'solver':'schism',
     'grid_file': GRID_FILE,
     'manning':.12,
     'windrot':0.00001,
     'tag':'schism',
     'start_date':'2018-10-1 0:0:0',
     'time_frame':'36H',
     'dem_source' : DEM_FILE,
     'meteo_source' : METEO_FILES_2,
     'engine':'cfgrib',
     'combine_forecast': True, #combine meteo
     'combine_by':'nested',
     'xr_kwargs': {'concat_dim':'step'},
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':1.5, 'nhot':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'nhot_write':108 }
    }



def schism(tmpdir):
    #initialize a model
    rpath = str(tmpdir)+'/schism/'
    case.update({'rpath':rpath+'20181001.00/'}) # use tmpdir for running the model

    b = pyPoseidon.model(**case)

    b.execute()

    # run the cast
    with open(rpath + '20181001.00/schism_model.json', 'rb') as f:
        info = pd.read_json(f,lines=True).T
        info[info.isnull().values] = None
        info = info.to_dict()[0]


    info.update({'path': rpath})

    #creating a time sequence of the runs
    start_date = pd.to_datetime('2018-10-1 0:0:0')
    end_date = pd.to_datetime('2018-10-2 0:0:0')
    date_list = pd.date_range(start_date,end_date, freq='12H')

    #append to dic
    info.update({'start_date':start_date,'end_date':end_date, 'dates' : date_list})

    #creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    #NOTE that the first folder is the fisrt run already perfomed!!
    folders = [datetime.datetime.strftime(x, '%Y%m%d.%H') for x in date_list]
    info.update({'folders':folders})

    #creating a sequence of folder from which we read the meteo.
    meteo = []
    for date in date_list:
        end_date= pd.to_datetime(date) + pd.to_timedelta(info['time_frame'])
        end_date = end_date.strftime(format='%Y-%m-%d %H:%M:%S')
        dr = pd.date_range(date, end_date, freq='12H')
        names = ['uvp_'+ datetime.datetime.strftime(x, '%Y%m%d%H') + '.grib' for x in dr]
        dur = [ (DATA_DIR / name).as_posix() for name in names ] 
        meteo.append(dur)

    info.update({'meteo_source':meteo})

    info['time_frame'] = len(folders)*[info['time_frame']]

    #set cast
    h = cast.cast(**info) # initialize

    h.run()

    # Run check case - Total duration
    check.update({'rpath':rpath+'check/'}) # use tmpdir for running the model

    c = pyPoseidon.model(**check)

    c.execute()

    # COMPARE
    folders = [info['path']+f for f in info['folders']]
    output = data.data(folders=folders,solver='schism')

    total = data.data(folders=[rpath+'check/'],solver='schism')


    rb = []
    for var in total.Dataset.data_vars:
        if not total.Dataset[var].equals(output.Dataset[var]):
            rb.append(var)

    print(rb)


#    flag = True TODO
#    for var in rb:
#        flag = False
#        mdif = np.abs(total.results.Dataset[var].values - output.results.Dataset[var].values).max()
#        if mdif < 1.e-14 :
#            flag = True
#    print(mdif)

    if (rb == ['zcor']) or rb==[]:
        return True
    else:
        return False

@pytest.mark.solvers
def test_answer(tmpdir):
    assert schism(tmpdir) == True
