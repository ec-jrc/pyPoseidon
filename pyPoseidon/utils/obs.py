"""
Observations management module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import requests, urllib
import datetime
from pydap.client import open_url
from dateutil.parser import parse
import numpy as np
import pandas as pd
import pkg_resources
import pyPoseidon
import os
import sys

#DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc')
DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    


class obs:

    def __init__(self,**kwargs):
        
        sdate = kwargs.get('start_date', None)
        edate = kwargs.get('end_date', None)
        sd = kwargs.get('sa_date', sdate)
        ed = kwargs.get('se_date', edate)
        
        self.sdate = pd.to_datetime(sd)
        self.edate = pd.to_datetime(ed)
        
        self.point = kwargs.get('point', None)
    
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)
        maxlat = kwargs.get('maxlat', None)
        
            
#        db = kwargs.get('filename', DATA_PATH+'ioc.csv')
        
#        ioc = pd.read_csv(db)
        critech = pd.read_csv(DATA_PATH+'critech.csv')
        
        critech.loc[:,['lon','lat']] = critech.loc[:,['lon','lat']].apply(pd.to_numeric)
     
        w = critech.loc[(critech['lon'] > minlon) & (critech['lon'] < maxlon) & (critech['lat'] > minlat) & (critech['lat'] < maxlat) \
                    & (pd.to_datetime(critech['Min. Time']) < self.sdate) & ~critech['Min. Time'].str.contains('Jan 0001') & ~critech['Min. Time'].str.contains('Jan 1970')]
        
        w.reset_index(inplace=True, drop=True)
        
        self.locations = w.copy() 
                                   

    def loc(self,name,**kwargs):
            
        point=self.locations[self.locations['Name'].str.contains(name)]['ID'].values[0]
        return self.webcritech(point=int(point))

    def iloc(self,idx,**kwargs):
            
        point=self.locations.iloc[idx,:]['ID']
        return self.webcritech(point=int(point))


    def webcritech(self,**kwargs):

        sdate = kwargs.get('start_date', self.sdate)
        edate = kwargs.get('end_date', self.edate)
        point = kwargs.get('point', self.point)

        pdate=min([self.edate+datetime.timedelta(hours=72),datetime.datetime.now()])        

        url='http://webcritech.jrc.ec.europa.eu/SeaLevelsDb/Home/ShowBuoyData?id={}&dateMin={}%2F{:02d}%2F{:02d}+{:02d}%3A{:02d}&dateMax={}%2F{:02d}%2F{:02d}+{:02d}%3A{:02d}&field=&options='\
                                 .format(point,sdate.year,sdate.month,sdate.day,sdate.hour,0,pdate.year,pdate.month,pdate.day,pdate.hour,0)
        
        #print(url)                                 
        response=requests.get(url)
        ls = response.content
        lp = ls.decode('utf-8').strip().split('\n')
        
  # get attrs
        try:
            #
            s1=lp[0].strip()
            at  = s1[s1.index('deviceID='):s1.index('(ID')]
            bname = at.strip().split('=')[1]
            at = s1[s1.index('(ID='):s1.index(')')]
            if len(at)==0: at = s1[s1.index('(ID='):]
            bid = at.strip().split('=')[1].strip(')')
            #
            s1=lp[1].strip()

            at = s1[s1.index('lat='):s1.index('lon')]
            lat = at.strip().split('=')[1]

            try:
                at = s1[s1.index('lon='):s1.index('hei')]
                lon = at.strip().split('=')[1]
            except:
                at = s1[s1.index('lon='):s1.index('Location')]
                lon = at.strip().split('=')[1]


            try:
                at = s1[s1.index('Location:'):s1.index('Extraction')]
            except:
                at = s1[s1.index('Location:'):s1.index('Last')]

            Location = at.strip().split(':')[1].strip('-').strip()
            Location = ''.join(c for c in Location if c not in '(){}<>\r\n#')


            try:
                at = s1[s1.index('Date:'):s1.index('Latency')]
                Last_Date = at.strip().split('Date:')[1].strip()
                at = s1[s1.index('Latency:'):]
                Latency = at.strip().split('Latency:')[1].strip()
                Extraction_date = np.nan
            except:
                at = s1[s1.index('Date:'):]
                Extraction_date = at.strip().split('Date:')[1].strip()
                Latency = np.nan
                Last_Date = np.nan

            #
            l = 2
            try:
                s1=lp[l].strip()
                at = s1[s1.index('Time:'):s1.index(' - Last')]
                Last_Time = at.strip().split('Time:')[1].strip()

                at = s1[s1.index('Value:'):]
                Last_Value = at.strip().split('Value:')[1].strip()
                Last_Value = ''.join(c for c in Last_Value if c not in '(){}<>\r\n#')

            except:
                l -= 1

            #
            l += 1
            s1=lp[l].strip()
            Server = s1.strip().split('=')[1]

            #
            l += 1
            s1=lp[l].strip()
            harmonics = pd.Series(s1.strip().split('=')[1].split('|'))
            harmonics = pd.DataFrame(harmonics.str.split(',', expand=True))
            harmonics.columns = ['h1','h2','h3']
            
            #
            l += 1
            s1 = lp[l]
            headers = ''.join(c for c in s1 if c not in '(){}<>\r\n#').split(',')
            headers = [x.strip() for x in headers]
            try:
                tg = pd.Series(lp[l+1:])
                tg = pd.DataFrame(tg.str.split(',', expand=True))
                tg.columns = headers
                tg[tg.columns[-1]] = tg[tg.columns[-1]].str.strip('\r')
                tg = tg.set_index(tg.columns[0])
                tg.index = pd.to_datetime(tg.index)
                tg = tg.apply(pd.to_numeric)
            except:
                tg = pd.DataFrame({'TimeUTC':np.nan, 'Level m':np.nan, 'Tide m':np.nan, 'Level-Tide m':np.nan}, index=[0])
            return tg
        except Exception as e:
            print(e)
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('problem with time series acquisition\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
            return None

        
#    def soest(self):
        
#        url = 'https://uhslc.soest.hawaii.edu/thredds/dodsC/uhslc/fdh/OS_UH-FDH329_20170628_D'
#        dataset = open_url(url)
        
#        t = dataset['time']
        
#        info = bunch(dataset.attributes['NC_GLOBAL'])
        
#        tref = t.attributes['units'].split()[-1]
        
#        self.tref = parse(tref)
        
        
#        time = [self.tref + datetime.timedelta(days = ta) for ta in t.data[:]]
        
        #find the index for the time frame we want
        #start_day =  
        
            
