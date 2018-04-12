#!/TEST/brey/miniconda3/envs/pyPoseidon/bin/python
# coding: utf-8
# ## hindcast example
from pyPoseidon.utils.cast import *

# ### setup
#read the info from the first run
with open('/COASTAL/DELFT3D/EUR/20100201.00/info.pkl', 'r') as f:
              info=pickle.load(f)
#define some info
info.update({'path':'/COASTAL/DELFT3D/EUR/', # The path of the project
     'case':'European 2010', # a reference tag
     })
#creating a time sequence of the runs
start_date = datetime.datetime(2010,2,1,0)
end_date = datetime.datetime(2010,3,1,12)
time_interval = datetime.timedelta(hours=12)
        
dt=(end_date-start_date).total_seconds()
ndt=dt/time_interval.total_seconds()
ndt=np.int(ndt)+1

date_list = [start_date + time_interval*x for x in range(ndt)]
#append to dic
info.update({'start_date':start_date,'end_date':end_date, 'time_interval':time_interval,'dates' : date_list})
#creating a sequence of folder to store the runs. In this case we name them after the date attribute.
#NOTE that the first folder is the fisrt run already perfomed!!
folders = [datetime.datetime.strftime(x, '%Y%m%d.%H') for x in date_list]
info.update({'folders':folders})

#creating a sequence of folder from which we read the meteo.
PATH='/COASTAL/meteo/'
meteo = [PATH+'{:04d}/{:02d}/{:02d}/'.format(x.year,x.month,x.day)+datetime.datetime.strftime(x, '%Y%m%d.%H')+'.tropical_cyclone.grib' for x in date_list]
info.update({'meteo_files':meteo})

# ### run the hindcast

#set cast
h = cast(**info) # initialize

h.run()
