import urllib2, urllib
import datetime
from pydap.client import open_url
from dateutil.parser import parse
import numpy as np
import pandas as pd

class obs:

    def __init__(self,**kwargs):
        
        self.sdate = kwargs.get('sdate', None)
        self.edate = kwargs.get('edate', None)
        self.point = kwargs.get('point', None)
        

    def getmes(self):

        pdate=min([self.edate+datetime.timedelta(hours=72),datetime.datetime.now()])


        url='http://webcritech.jrc.ec.europa.eu/SeaLevelsDb/Home/ShowBuoyData?id={}&dateMin={}%2F{:02d}%2F{:02d}+{:02d}%3A{:02d}&dateMax={}%2F{:02d}%2F{:02d}+{:02d}%3A{:02d}&field=&options='\
                                 .format(self.point,self.sdate.year,self.sdate.month,self.sdate.day,self.sdate.hour,0,pdate.year,pdate.month,pdate.day,pdate.hour,0)
        response=urllib2.urlopen(url)
        ls=response.readlines()
        lp=[elem.strip().split(',')  for elem in ls]
  # get name id
        try:
            lp0=''.join(lp[0])
            bname=lp0.split('ID=')[1].strip('(=)').strip()
            bid=lp0.split('ID=')[2].strip('(=)').strip()
        except:
            lp1=''.join(lp[1])
            bname=lp1.split('ID=')[1].strip('(=)').strip()
            bid=lp1.split('ID=')[2].strip('(=)').strip()

  # get lat lon
        c=[a.split(' ') for a in lp[1]][0]
        if 'lat=' in c[2]: 
           plat=c[2].split('=')[1]
           idt=6
        else:
           c=[a.split(' ') for a in lp[2]][0]
           if 'lat=' in c[2]: plat=c[2].split('=')[1]
           idt=7

        if 'lon=' in c[3]: plon=c[3].split('=')[1]

        rt=[]
        vt=[]
        for a,b,c,d in lp[idt:]:
            rt.append(datetime.datetime.strptime(a,'%d %b %Y %H:%M:%S'))
            vt.append([b,c,d])

        tg = pd.DataFrame(np.column_stack([rt,vt]),columns = ['time', 'Total Sea Level', 'Tide', 'Storm Surge'])
        tg = tg.set_index(['time'])
        tg = tg.apply(pd.to_numeric)
        
        self.data = tg

        
    def soest(self):
        
        url = 'https://uhslc.soest.hawaii.edu/thredds/dodsC/uhslc/fdh/OS_UH-FDH329_20170628_D'
        dataset = open_url(url)
        
        t = dataset['time']
        
        info = bunch(dataset.attributes['NC_GLOBAL'])
        
        tref = t.attributes['units'].split()[-1]
        
        self.tref = parse(tref)
        
        
        time = [self.tref + datetime.timedelta(days = ta) for ta in t.data[:]]
        
        #find the index for the time frame we want
        #start_day =  
        
            
