#!/TEST/brey/miniconda3/envs/pyPoseidon/bin/python
# coding: utf-8
# ## An example workflow
from pyPoseidon.model import *
from pyPoseidon.dem import *
from pyPoseidon.meteo import *
import datetime

#define in a dictionary the properties of the model..
dic={'minlon':-35., # lat/lon window
     'maxlon':42.,
     'minlat':25.05,
     'maxlat':76.5,
     'solver':'d3d',
     'resolution':0.1, #grid resoltuion 
     'step':60, # step for output of map field in d3d
     'rstep':60 * 12, # step for output of restart file    
     'start_date':'2010-2-1',
     'time_frame':'12H',
     'meteo':'ecmwf_oper',
     'dem': 'emodnet',
     'dpath' : '/COASTAL/BATHYMETRY/EMODNET/TEMP_bathymetry.nc',
     'exec':'/TEST/brey/DELFT3D/7545/bin/lnx64/', #exec folder of solver 
     'ncores': 32, #number of cores
     'rpath':'/COASTAL/DELFT3D/EUR/20100201.00/', #location of calc folder
    }

# #### Local operational ECMWF files

if 'time_frame' in dic.keys(): end_date= pd.to_datetime(dic['start_date']) + pd.to_timedelta(dic['time_frame'])
dic.update({'end_date':end_date.strftime(format='%Y-%m-%d')})

dr = pd.date_range(dic['start_date'],dic['end_date'], freq='12H')

#creating a sequence of folder from which we read the meteo.
PATH='/COASTAL/meteo/'#Path to meteo files
folders = [datetime.datetime.strftime(x, '%Y%m%d.%H') for x in dr]
meteo = [PATH+'{:04d}/{:02d}/{:02d}/'.format(x.year,x.month,x.day)+datetime.datetime.strftime(x, '%Y%m%d.%H')+'.tropical_cyclone.grib' for x in dr]

# #### update dictionary

dic.update({'mpaths':meteo,'ft1':0,'ft2':12})


# ## Initialize

#initialize a model
b = model(**dic)

# ### set it up
b.set(**dic) #set it up 

# #### Optional adjust the wet area based on a coastline shapefile

#b.impl.dem.impl.adjust('/COASTAL/COASTLINES/naturalearth/ne_50m_coastline.shp')#,nc=20) 

# ## Save to folder for execution 

#set the run by saving the files
b.output()

# save model info for further use
b.save()

# ### execute

#execute interactively
b.run()

#execute batch
#b.run()
