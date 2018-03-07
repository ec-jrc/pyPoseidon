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
import xarray as xr
import pandas as pd
import holoviews as hv
from cartopy import crs
import geoviews as gv
import importlib
from pyPoseidon.utils.get_value import get_value
from dateutil import parser


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

            date=grib_get(gid, 'date')
            dataTime=grib_get(gid, 'dataTime')
            stepRange=grib_get(gid, 'stepRange')

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

            return date,dataTime,stepRange,name,dat,lon,lat


class meteo:
    impl=None
    def __init__(self,**kwargs):
        msource = kwargs.get('meteo', None)
        if msource == 'ecmwf_oper' :
            self.impl = ecmwf_oper(**kwargs)
        else:
            self.impl = gfs(**kwargs)
                
    
class ecmwf_oper(meteo):   
        
    def __init__(self,**kwargs):
    
      filenames = kwargs.get('mpaths', {})
      ft1 = kwargs.get('ft1', None)
      ft2 = kwargs.get('ft2', None)
      dft = kwargs.get('dft', 1)
      
      ft2 = ft2 + 1 # for range
                      
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
    
      nt1=3*ft1
      nt2=3*ft2
      step = 3*dft            
      
      # read grib file and append to xarray

      pt=[]
      ut=[]
      vt=[]
      tt=[]

      for filename in filenames:

     	try: 
       		f = open(filename)
      	except:
        	print 'no file {}'.format(filename)
         	sys.exit(1)


    	for it in xrange(nt1):
            gid = grib_new_from_file(f)#,headers_only = True)
            grib_release(gid)

      #--------------------------------------------------------------------- 
      	sys.stdout.flush()
      	sys.stdout.write('\n')
      	sys.stdout.write('extracting meteo from {}\n'.format(filename))
      	sys.stdout.flush()
      #---------------------------------------------------------------------      

      	indx = [[i,i+1,i+2] for i in np.arange(nt1,nt2,step)]
      	flag = [item for sublist in indx for item in sublist]


      	try:
          for it in tqdm(range(nt1,nt2)): 
            
            if it not in flag : # skip unwanted timesteps
                gid = grib_new_from_file(f)#,headers_only = True)
                grib_release(gid)
                continue           
            
            date,dataTime,stepRange,name,varin,ilon,ilat=getd(f,it)    

            timestamp =  datetime.datetime.strptime(str(date),'%Y%m%d')+datetime.timedelta(hours=dataTime/100) 

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
                     tt.append(timestamp+datetime.timedelta(hours=int(stepRange)))
            elif name == '10u':
                     ut.append(data)
            elif name == '10v':
                     vt.append(data)


    # END OF FOR

        except Exception as e:
          print e
          print 'ERROR in meteo file {}'.format(date)

        f.close()

      met = xr.Dataset({'msl': (['time', 'latitude', 'longitude'],  np.array(pt)), 
                          'u10': (['time', 'latitude', 'longitude'], np.array(ut)),   
                          'v10': (['time', 'latitude', 'longitude'], np.array(vt)),   
                          'lons': (['x', 'y'], lons),   
                          'lats': (['x', 'y'], lats)},   
                          coords={'longitude': ('longitude', lons[0,:]),   
                                  'latitude': ('latitude', lats[:,0]),   
                                  'time': tt })   
#                        'time': pd.date_range(date+datetime.timedelta(hours=ft1), periods=ft2-ft1, freq='{}H'.format(dft))})   
#                        'reference_time': date })
    
      self.p = np.array(pt)
      self.u = np.array(ut)
      self.v = np.array(vt)
      self.uvp = met
      self.lats = lats
      self.lons = lons   
      self.ft1 = ft1
      self.ft2 = ft2 
      self.dft = dft   
      
      self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

      self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('meteo done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
    
        
            
      

class gfs(meteo):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
      ts = kwargs.get('start_date', None)
      te = kwargs.get('end_date', None)
      
      ts = pd.to_datetime(ts)
      te = pd.to_datetime(te)     
         
      url = kwargs.get('url', 'https://bluehub.jrc.ec.europa.eu/erddap/griddap/NCEP_Global_Best')
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting meteo from {}\n'.format(url))
      sys.stdout.flush()
      #---------------------------------------------------------------------      

      data = xr.open_dataset(url)    
      
      if minlon < data.geospatial_lon_min : minlon = minlon + 360.
      
      if maxlon < data.geospatial_lon_min : maxlon = maxlon + 360.
      
      if minlon > data.geospatial_lon_max : minlon = minlon - 360.
      
      if maxlon > data.geospatial_lon_max : maxlon = maxlon - 360.
            
      if ts < parser.parse(data.attrs['time_coverage_start']).replace(tzinfo=None) :
          sys.stdout.flush()
          sys.stdout.write('\n')
          sys.stdout.write('time frame not available\n')
          sys.stdout.write('coverage between {} and {} \n'.format(data.attrs['time_coverage_start'],data.attrs['time_coverage_end']))
          sys.stdout.flush()
          sys.exit(1)
          
      if te > parser.parse(data.attrs['time_coverage_end']).replace(tzinfo=None) :
          sys.stdout.flush()
          sys.stdout.write('\n')
          sys.stdout.write('time frame not available\n')
          sys.stdout.write('coverage between {} and {} \n'.format(data.attrs['time_coverage_start'],data.attrs['time_coverage_end']))
          sys.stdout.flush()
          sys.exit(1)
      
      tslice=slice(ts, te)
    
      i0=np.abs(data.longitude.data-minlon).argmin()
      i1=np.abs(data.longitude.data-maxlon).argmin()

      
      j0=np.abs(data.latitude.data-minlat).argmin()
      j1=np.abs(data.latitude.data-maxlat).argmin()

      if i0 > i1 :

          sh = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(longitude=slice(i0,data.longitude.size),latitude=slice(j0,j1))
              .sel(time=tslice)
              )
          sh.longitude.values = sh.longitude.values -360.

          sh1 = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(longitude=slice(0,i1),latitude=slice(j0,j1))
              .sel(time=tslice)
              )
              
          tot = xr.concat([sh,sh1],dim='longitude')
          
      else:            

          tot = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(longitude=slice(i0,i1),latitude=slice(j0,j1))
              .sel(time=tslice)
              )
              
      if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
          c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
          tot['longitude'] = tot['longitude'] + c*360.
      
      xx,yy=np.meshgrid(tot.longitude,tot.latitude)
      
      other = xr.Dataset({'lons': (['x', 'y'], xx), 
                          'lats': (['x', 'y'], yy) })
      
      self.uvp = xr.merge([tot, other])
      
      self.uvp.rename({'prmslmsl':'msl','ugrd10m':'u10','vgrd10m':'v10'}, inplace=True)
      
      self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

      self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('meteo done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
        
        
class generic(meteo):
    
    def __init__(self,**kwargs):
    
        filename = kwargs.get('filename', None)
        
        self.uvp = xr.open_dataset(filename)
         
        self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

        self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
                
    
    def output(self,solver=None,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
