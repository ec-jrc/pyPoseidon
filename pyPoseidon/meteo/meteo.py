import numpy as np
import datetime
import glob
import sys
from time import sleep
from copy import deepcopy
import os
from gribapi import *
from redtoreg import _redtoreg
from pygrib import gaulats
from tqdm import tqdm
import time


def gridd(lon1,lat1,lon2,lat2,nlats):

            #   lon1, lat1 = self.longitude_first_gridpoint, self.latitude_first_gridpoint
            #   lon2, lat2 = self.longitude_last_gridpoint, self.latitude_last_gridpoint
            #   nlats = self.points_in_y_direction
                # ECMWF 'reduced' gaussian grid.
                nlons = 2*nlats
                delon = 360./nlons
            #   lons = np.arange(lon1,lon2,delon)
                lons = np.linspace(lon1,lon2,nlons)
                # compute gaussian lats (north to south)
                lats = gaulats(nlats)
                if lat1 > lat2 :
                   lats = lats[::-1]
              # lons = lons[::-1]
                lons,lats = np.meshgrid(lons,lats) # make 2-d arrays

                return lons,lats


def getd(f,t):
            gid = grib_new_from_file(f)#,headers_only = True)
            if gid is None: 
                print 'time = {}, gid = None'.format(t)
                sys.exit(1)

            name=grib_get(gid, 'shortName')
            mv=grib_get(gid,'missingValue')

            lonfgp=grib_get(gid,'longitudeOfFirstGridPointInDegrees')
            latfgp=grib_get(gid,'latitudeOfFirstGridPointInDegrees')
            lonlgp=grib_get(gid,'longitudeOfLastGridPointInDegrees')
            latlgp=grib_get(gid,'latitudeOfLastGridPointInDegrees')

            if grib_get(gid,'gridType') == 'regular_gg':

              Ni=grib_get(gid,'Ni')
              Nj=grib_get(gid,'Nj')
              lat=grib_get_array(gid,'latitudes')
              lat=lat.reshape(Nj,Ni)
              lat=np.flipud(lat)
              lon=grib_get_array(gid,'longitudes')
              lon=lon.reshape(Nj,Ni)

              values=grib_get_values(gid)
              dat=values.reshape(Nj,Ni)
              dat=np.flipud(dat)
          
            elif grib_get(gid,'gridType') == 'reduced_gg' :

              ss=grib_get_array(gid,'pl')  # lons per lat for the reduced_gg grid
              lon,lat = gridd(lonfgp,latfgp,lonlgp,latlgp,ss.size)

              values=grib_get_values(gid)
              ny=2*np.size(ss)

              dat=_redtoreg(ny,ss,values,mv)
              dat=np.flipud(dat)

            grib_release(gid)

            return name,dat,lon,lat


class meteo:
    impl=None
    def __init__(self,**kwargs):
        msource = kwargs.get('meteo', None)
        if msource == 'ecmwf' :
            self.impl = ecmwf(**kwargs)
                
    
class ecmwf(meteo):   
        
    def __init__(self,**kwargs):
    
      filename = kwargs.get('mpath', {})
      ft1 = kwargs.get('ft1', None)
      ft2 = kwargs.get('ft2', None)
      dft = kwargs.get('dft', None)
      
      ft2 = ft2 + 1 # for range
                      
      minlon = kwargs.get('lon0', None)
      maxlon = kwargs.get('lon1', None)
      minlat = kwargs.get('lat0', None)
      maxlat = kwargs.get('lat1', None) 
    
      nt1=3*ft1
      nt2=3*ft2
      step = 3*dft            
      
      # read grib file

      try: 
       f = open(filename)
      except:
         print 'no file {}'.format(filename)
         sys.exit()

      for it in xrange(nt1):
            gid = grib_new_from_file(f)#,headers_only = True)
            grib_release(gid)

      pt=[]
      ut=[]
      vt=[]
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting meteo\n')
      sys.stdout.flush()
      #---------------------------------------------------------------------      
      

      mxv=nt2-nt1-1
      try:
        for it in tqdm(range(nt1,nt2,step)): 
            
            name,varin,ilon,ilat=getd(f,it)        

            lon=ilon[0,:]
            lat=ilat[:,0]

        # get sea level pressure and 10-m wind data.
        
        # shift grid according to minlon
            if minlon < 0. :
               lon=lon-180.
               zlon=lon.shape[0]
               varin_ = np.hstack([varin[:,zlon/2:],varin[:,0:zlon/2]])
               varin  = varin_

            i1=np.abs(lon-minlon).argmin()-2
            i2=np.abs(lon-maxlon).argmin()+2
            j1=np.abs(lat-minlat).argmin()-2
            j2=np.abs(lat-maxlat).argmin()+2
    
            if i1 < 0 : i1 = 0 # fix limits
            if i2 > lon.shape[0] : i2 = lon.shape[0]   
            if j1 < 0 : j1 = 0
            if j2 > lat.shape[0]: j2 = lat.shape[0]
            

            lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
            data = deepcopy(varin[j1:j2,i1:i2])

        # mask the window

            if name == 'msl' : 
                     pt.append(data)
            elif name == '10u':
                     ut.append(data)
            elif name == '10v':
                     vt.append(data)


    # END OF FOR
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('meteo done\n')
        sys.stdout.flush()
        #--------------------------------------------------------------------- 

      except Exception as e:
        print e
        print 'ERROR in meteo input'

      f.close()
    
      self.p = np.array(pt)
      self.u = np.array(ut)
      self.v = np.array(vt)
      self.lats = lats
      self.lons = lons   
      self.ft1 = ft1
      self.ft2 = ft2    
      
      
      
