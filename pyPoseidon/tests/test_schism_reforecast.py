import pytest
import pyPoseidon
from pyPoseidon.utils import cast, data
import pyPoseidon.meteo as pm
import json
import pandas as pd
import datetime
import os
import numpy as np
import xarray as xr

from . import DATA_DIR


GRID_FILE = (DATA_DIR / 'hgrid.gr3').as_posix()
DEM_FILE = (DATA_DIR / 'dem.nc').as_posix()
METEO_FILES_1 = [(DATA_DIR / 'uvp_2018100100.grib').as_posix()]
METEO_FILES_2 = [(DATA_DIR / name).as_posix() for name in ('uvp_2018100100.grib', 'uvp_2018100112.grib', 'uvp_2018100200.grib', 'uvp_2018100212.grib')]


#define in a dictionary the properties of the model..
case={'solver':'schism',
     'grid_file': GRID_FILE,
     'manning':.12,
     'windrot':0.00001,
     'tag':'schism',
     'start_date':'2018-10-1 0:0:0',
     'time_frame':'24H',
     'dem_source' : DEM_FILE,
     'engine':'passthrough',
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':1., 'nhot':1, 'ihot':0,'nspool':9, 'ihfskip':36, 'nhot_write':108 }
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
     'engine':'passthrough',
     'ncores': 4 , #number of cores
     'update':['all'], #update only meteo, keep dem
     'parameters':{'dt':400, 'rnday':1.5, 'nhot':0, 'ihot':0,'nspool':9, 'ihfskip':36, 'nhot_write':108 }
    }



def schism(tmpdir):
    #initialize a model
    rpath = str(tmpdir)+'/schism/'
    case.update({'rpath':rpath+'20181001.00/'}) # use tmpdir for running the model

    #creating a time sequence of the runs
    start_date = pd.to_datetime('2018-10-1 0:0:0')
    end_date = pd.to_datetime('2018-10-2 0:0:0')
    date_list = pd.date_range(start_date,end_date, freq='12H')


    m0 = pm.meteo(meteo_source=METEO_FILES_1,engine='cfgrib')
    
    case.update({'meteo_source':m0.Dataset}) 

    b = pyPoseidon.model(**case)

    b.execute()

    # run the cast
    with open(rpath + '20181001.00/schism_model.json', 'rb') as f:
        info = pd.read_json(f,lines=True).T
        info[info.isnull().values] = None
        info = info.to_dict()[0]


    info.update({'path': rpath})


    #append to dic
    info.update({'start_date':start_date,'end_date':end_date, 'dates' : date_list})

    #creating a sequence of folder to store the runs. In this case we name them after the date attribute.
    #NOTE that the first folder is the fisrt run already perfomed!!
    folders = [datetime.datetime.strftime(x, '%Y%m%d.%H') for x in date_list]
    info.update({'folders':folders})

    #creating a sequence of folder from which we read the meteo.
    meteo = [m0.Dataset]
    for date in date_list[1:]:
        end_date= pd.to_datetime(date) + pd.to_timedelta(info['time_frame'])
        end_date = end_date.strftime(format='%Y-%m-%d %H:%M:%S')
        dr = [date - pd.to_timedelta('12H'), date]
        names = ['uvp_'+ datetime.datetime.strftime(x, '%Y%m%d%H') + '.grib' for x in dr]
        dur = [ (DATA_DIR / name).as_posix() for name in names] 
        m1 = pm.meteo(meteo_source=dur[0],engine='cfgrib')
        m2 = pm.meteo(meteo_source=dur[1],engine='cfgrib')
        w1 = m1.Dataset.isel(time=slice(12,13))
        w2 = m2.Dataset.isel(time=slice(1,None)) # note that we keep the 12 hour from the previous file
        mf = xr.combine_by_coords([w1,w2])
        meteo.append(mf)

    info.update({'meteo_source':meteo})

    info['time_frame'] = len(folders)*[info['time_frame']]

    #set cast
    h = cast.cast(**info) # initialize

    h.run()

    # Run check case - Total duration
    check.update({'rpath':rpath+'check/'}) # use tmpdir for running the model

    # Combine meteo appropriately

    m1 = pm.meteo(meteo_source=METEO_FILES_2[0],engine='cfgrib')
    m2 = pm.meteo(meteo_source=METEO_FILES_2[1],engine='cfgrib')
    m3 = pm.meteo(meteo_source=METEO_FILES_2[2],engine='cfgrib')
    m4 = pm.meteo(meteo_source=METEO_FILES_2[3],engine='cfgrib')

    # extract correct chunk

    w1 = m1.Dataset.isel(time=slice(0,13))
    w2 = m2.Dataset.isel(time=slice(1,13)) # note that we keep the 12 hour from the previous file
    w3 = m3.Dataset.isel(time=slice(1,13))
    w4 = m4.Dataset.isel(time=slice(1,13))

    #combine
    meteo = xr.combine_by_coords([w1,w2,w3,w4])
    #saving
    check.update({'meteo_source' : meteo})
    
    c = pyPoseidon.model(**check)

    c.execute()

    # COMPARE
    folders = [info['path']+f for f in info['folders']]
    output = data.data(folders=folders,solver='schism')

    total = data.data(folders=[rpath+'check/'],solver='schism')

    r = output.Dataset.isel(time=slice(0,36))
    

    rb = []
    for var in total.Dataset.data_vars:
        if not total.Dataset[var].equals(r[var]):
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

@pytest.mark.schism
def test_answer(tmpdir):
    assert schism(tmpdir) == True
