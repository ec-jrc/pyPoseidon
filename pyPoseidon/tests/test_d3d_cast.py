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
GRIB_FILES_2 = [(DATA_DIR / filename).as_posix() for filename in ('uvp_2018100100.grib', 'uvp_2018100112.grib', 'uvp_2018100200.grib', 'uvp_2018100212.grib')]
DEM_FILE = (DATA_DIR / "dem.nc").as_posix()

#define in a dictionary the properties of the model..
case={'geometry':{'lon_min' : -30,'lon_max' : -10.,'lat_min' : 60.,'lat_max' : 70.},
     'start_date':'2018-10-1',
     'time_frame':'12H',
     'solver':'d3d',
     'resolution':0.1, #grid resoltuion
     'map_step':60, # step for output of map field in d3d
     'restart_step':720, # when to output restart file
     'ncores': NCORES , #number of cores
     'meteo_source' : GRIB_FILES_1,
     'meteo_engine':'cfgrib',
     'meteo_merge': 'last', #combine meteo
     'meteo_combine_by':'nested',
     'meteo_xr_kwargs': {'concat_dim':'step'},
     'dem_source' : DEM_FILE,
#     'update':['all'] # optional to select update quantities
    }

check={'geometry':{'lon_min' : -30,'lon_max' : -10.,'lat_min' : 60.,'lat_max' : 70.},
     'start_date':'2018-10-1',
     'time_frame':'36H',
     'solver':'d3d',
     'resolution':0.1, #grid resoltuion
     'map_step':60, # step for output of map field in d3d
     'restart_step':720, # when to output restart file
     'ncores': NCORES , #number of cores
     'dem_source' : DEM_FILE,
     'meteo_source' : GRIB_FILES_2,
     'meteo_engine':'cfgrib',
     'meteo_merge': 'last', #combine meteo
     'meteo_combine_by':'nested',
     'meteo_xr_kwargs': {'concat_dim':'step'},
#     'update':['all'] # optional to select update quantities
    }


def d3d(tmpdir):
    #initialize a model
    rpath = str(tmpdir) + '/d3d/'
    case.update({'rpath':rpath + '/20181001.00/'}) # use tmpdir for running the model
    b = pyPoseidon.model(**case)

    b.execute()

    #creating a time sequence of the runs
    start_date = pd.to_datetime('2018-10-1 0:0:0')
    end_date = pd.to_datetime('2018-10-2 0:0:0')
    date_list = pd.date_range(start_date,end_date, freq='12H')


    #creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    #NOTE that the first folder is the fisrt run already perfomed!!
    rpaths = [rpath + datetime.datetime.strftime(x, '%Y%m%d.%H') +'/' for x in date_list]

    #set meteo files
    meteo = []
    for date in date_list:
        end_date= pd.to_datetime(date) + pd.to_timedelta('12H')
        end_date = end_date.strftime(format='%Y-%m-%d %H:%M:%S')
        dr = pd.date_range(date, end_date, freq='12H')
        names = ['uvp_'+ datetime.datetime.strftime(x, '%Y%m%d%H') + '.grib' for x in dr]
        dur = [ (DATA_DIR / name).as_posix() for name in names ]
        meteo.append(dur)
        
        
    #set cast
    for l in range(len(rpaths)-1):
        h = cast.cast(solver='d3d',model=b,ppath=rpaths[l],cpath=rpaths[l+1],meteo=meteo[l+1], date=date_list[l+1])
        h.set(execute=True) # execute
    
    # Run check case - Total duration
    check.update({'rpath':rpath+'check/'}) # use tmpdir for running the model
    c = pyPoseidon.model(**check)
    c.execute()


    # COMPARE
    output = data.data(folders=rpaths,solver='d3d')
    total = data.data(folders=[rpath+'check/'],solver='d3d')

    test = True
    rb=[]
    for var in total.Dataset.data_vars:
        if not total.Dataset[var].equals(output.Dataset[var]):
            rb.append(var)
            if np.abs(total.Dataset[var].values-output.Dataset[var].values).max() > 1.e-6 : test = False

    print(rb)
    return test


#@pytest.mark.delft
def test_answer(tmpdir):
    assert d3d(tmpdir) == True
