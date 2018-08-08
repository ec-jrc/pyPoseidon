#!/Users/brey/miniconda3/envs/dev/bin/python
# ## hindcast example
from pyPoseidon.utils import cast, data, point
import pickle
import pandas as pd
import datetime
# ### setup
#read the info from the first run
with open('/Users/brey/Downloads/EUR/D3D/ERAI/201701/gebco_eur_01_info.pkl', 'r') as f:
              info=pickle.load(f)

#define some info
info.update({'path':'/Users/brey/Downloads/EUR/D3D/ERAI/', # The path of the project
     'case':'European 2010', # a reference tag
     })

start_date = info['start_date']
end_date = pd.to_datetime('2017-5-1 0:0:0')

s = pd.Series(index=pd.date_range(start_date, end_date))
df = s.resample('MS').size().to_period('m').rename_axis('Month').reset_index(name='NumDays')

#creating a time sequence of the runs
date_list = pd.date_range(start_date,end_date, freq='1M') - pd.offsets.MonthBegin(1, normalize=True)

dates = [pd.to_datetime(x) for x in date_list]

#creating a sequence of folder to store the runs. In this case we name them after the date attribute.
#NOTE that the first folder is the fisrt run already perfomed!!
folders = [x.strftime('%Y%m') for x in date_list]

time_frames = [pd.to_timedelta('{}D'.format(x)) for x in df.NumDays.values[:-1]]

#append to dic
info.update({'start_date':start_date,'end_date':end_date, 'dates' : dates, 'time_frames':time_frames, 'folders':folders})

#creating a sequence of folder from which we read the meteo.
meteo = []
PATH='/Users/brey/DATA/ERAI/'
for date in date_list:
    year = int(pd.to_datetime(date).strftime('%Y'))
    dur = [PATH+'eraInterim_{:04d}.grib'.format(year)]
    meteo.append(dur)

info.update({'meteo_files':meteo})

# ### run the hindcast

#set cast
h = cast.cast(**info) # initialize

h.run()


