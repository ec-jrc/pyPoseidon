import pytest
import pyPoseidon
from pyPoseidon.utils import cast, data
import json
import pandas as pd
import datetime
import os
import numpy as np
import multiprocessing

NCORES = max(1, multiprocessing.cpu_count() - 1)

from . import DATA_DIR

GRIB_FILES_1 = [(DATA_DIR / filename).as_posix() for filename in ('uvp_2018100100.grib', 'uvp_2018100112.grib')]
GRIB_FILES_2 = [(DATA_DIR / filename).as_posix() for filename in ('uvp_2018100100.grib', 'uvp_2018100112.grib', 'uvp_2018100200.grib')]
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

#define in a dictionary the properties of the model..
case1={'geometry':{'lon_min' : -30,'lon_max' : -10.,'lat_min' : 60.,'lat_max' : 70.},
     'start_date':'2018-10-1',
     'time_frame':'12H',
     'solver':'d3d',
     'resolution':0.1, #grid resoltuion
     'map_step':60, # step for output of map field in d3d
     'restart_step':720, # when to output restart file
     'ncores': NCORES , #number of cores
     'meteo_source' : GRIB_FILES_1,
     'engine':'cfgrib',
     'combine_forecast' : True,
     'combine_by':'nested',
     'xr_kwargs': {'concat_dim':'step'},
     'dem_source' : DEM_FILE,
#     'update':['all'] # optional to select update quantities
    }

case2={'geometry':{'lon_min' : -30,'lon_max' : -10.,'lat_min' : 60.,'lat_max' : 70.},
     'start_date':'2018-10-1',
     'time_frame':'24H',
     'solver':'d3d',
     'resolution':0.1, #grid resoltuion
     'map_step':60, # step for output of map field in d3d
     'restart_step':720, # when to output restart file
     'ncores': NCORES , #number of cores
     'dem_source' : DEM_FILE,
     'meteo_source' : GRIB_FILES_2,
     'engine':'cfgrib',
     'combine_forecast' : True,
     'combine_by':'nested',
     'xr_kwargs': {'concat_dim':'step'},
#     'update':['all'] # optional to select update quantities
    }


def d3d(tmpdir,dic):
    #initialize a model
    rpath = str(tmpdir)
    dic.update({'rpath':rpath + '/20181001.00/'}) # use tmpdir for running the model
    b = pyPoseidon.model(**dic)

    b.execute()
    # Cast
    #read the info from the first run
    with open(rpath+'/20181001.00/d3d_model.json', 'rb') as f:
        info = pd.read_json(f,lines=True).T
        info[info.isnull().values] = None
        info = info.to_dict()[0]

    info.update({'path':rpath}) # The path of the project

    #creating a time sequence of the runs
    start_date = pd.to_datetime('2018-10-1 0:0:0')
    end_date = pd.to_datetime('2018-10-1 12:0:0')
    date_list = pd.date_range(start_date,end_date, freq='12H')
    #append to dic
    info.update({'start_date':start_date,'end_date':end_date, 'dates' : date_list})

    #creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    #NOTE that the first folder is the fisrt run already perfomed!!
    folders = [datetime.datetime.strftime(x, '%Y%m%d.%H') for x in date_list]
    info.update({'folders':folders})

    #set meteo files
    meteo = []
    for date in date_list:
        end_date= pd.to_datetime(date) + pd.to_timedelta(info['time_frame'])
        end_date = end_date.strftime(format='%Y-%m-%d %H:%M:%S')
        dr = pd.date_range(date, end_date, freq='12H')
        dur = [(DATA_DIR / ('uvp_' + datetime.datetime.strftime(x, '%Y%m%d%H') + '.grib')).as_posix() for x in dr]
        meteo.append(dur)
    info.update({'meteo_source':meteo})
    print(meteo)

    info.update({'time_frame' : len(folders)*[info['time_frame']]})


    h = cast.cast(**info) # initialize
    h.run()
    # combine output
    folders = [info['path']+'/'+f for f in info['folders']]
    res = data.data(folders=folders,solver='d3d')


    # check single run
    case2.update({'rpath':rpath + '/combined/'})
    a = pyPoseidon.model(**case2)
    a.execute()
    out = data.data(**case2)

    test = True
    for var in out.Dataset.data_vars:
        if not out.Dataset[var].equals(res.Dataset[var]):
            if np.abs(out.Dataset[var].values-res.Dataset[var].values).max() > 1.e-6 : test = False


    return test


@pytest.mark.delft
@pytest.mark.parametrize('case', [case1])
def test_answer(tmpdir, case):
    assert d3d(tmpdir,case) == True
