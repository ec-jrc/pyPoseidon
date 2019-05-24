# -*- coding: utf-8 -*- 
"""
Meteo module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

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
import time
import xarray as xr
import pandas as pd
import holoviews as hv
from cartopy import crs
import geoviews as gv
import importlib
from pyPoseidon.utils.get_value import get_value
from dateutil import parser
from pyresample import bilinear, geometry, kd_tree
from pyresample import utils
import pyproj


def gridd(lon1,lat1,lon2,lat2,nlats):

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


def getd(gid):
           # gid = grib_new_from_file(f)#,headers_only = True)

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

            elif grib_get(gid,'gridType') == 'regular_ll':

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
        
            grib_release(gid)

            return name,dat,lon,lat


class meteo:
    impl=None
    def __init__(self,**kwargs):
        msource = kwargs.get('meteo', None)
        if msource == 'ecmwf_oper' :
            self.impl = grib_cfgrib(**kwargs)
        elif msource == 'hnms_oper' :
            self.impl = grib_pynio_HNMS(**kwargs)
        elif msource == 'am_oper' :
            self.impl = grib_pynio_AM(**kwargs)
        elif msource == 'generic' :
            self.impl = generic(**kwargs)
        elif msource == 'gfs_oper' :
            self.impl = gfs_oper(**kwargs)
        elif msource in ['erai','era5'] :
            self.impl = grib_cfgrib(**kwargs)
            
        
        else:
            self.impl = gfs_erdap(**kwargs)


class combine_meteo(meteo):
    
    def __init__(self,**kwargs):
        
        filenames = kwargs.get('mpaths', {})
                
        dats = []
        for filename in filenames[:-1]:
            
            kwargs['mpaths'] = filename
            ds = grib_cfgrib(**kwargs)
            dats.append(ds.uvp)
        
        kwargs['ft2'] += 1
        kwargs['mpaths'] = filenames[-1]
        ds = grib_cfgrib(**kwargs)
        dats.append(ds.uvp)
        
        self.uvp = xr.concat(dats,dim='time')
        
        
    def output(self,solver=None,**kwargs):
    
        path = get_value(self,kwargs,'rpath','./') 
           
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
           
        s = getattr(model,solver) # get solver class
           
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
   
                      

class grib_cfgrib(meteo):
    
    def __init__(self,**kwargs):
    
    
        filenames = kwargs.get('mpaths', {})
        ft1 = kwargs.get('ft1', None)
        ft2 = kwargs.get('ft2', None)
        dft = kwargs.get('dft', None)
        engine = kwargs.get('engine', 'cfgrib')
        backend_kwargs = kwargs.get('backend_kwargs', {'indexpath':''})
    
        start_date = kwargs.get('start_date', None)
        try:
            start_date = pd.to_datetime(start_date)
        except:
            pass
    
        if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
        else:
          end_date = kwargs.get('end_date', None)
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass
        
        ts = pd.to_datetime(start_date)
        te = pd.to_datetime(end_date)     
        
                          
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
      
      #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo')
        sys.stdout.flush()
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, backend_kwargs=backend_kwargs)    
        
        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
        
        time_coord = data.msl.dims[0]
        
        if time_coord != 'time':
            if 'time' in data.coords.keys():
                data = data.rename({'time':'rtime'})
                data = data.rename({time_coord:'time'})
                data = data.assign_coords(time=data.valid_time)
                
        if not minlon : minlon = data.longitude.data.min()
        if not maxlon : maxlon = data.longitude.data.max()
        if not minlat : minlat = data.latitude.data.min()
        if not maxlat : maxlat = data.latitude.data.max()
        
        
        if minlon < data.longitude.data.min() : minlon = minlon + 360.
      
        if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.
      
        if minlon > data.longitude.data.max() : minlon = minlon - 360.
      
        if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
                
        
        try:
            data = data.assign_coords(time=data.valid_time)
        except:
            sys.exit(1)
            
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        if ft1 :
            ts = data.time.data[ft1]
        
        if ft2 :
            te = data.time.data[ft2]
        
        
        if ts < data.time.min() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
          
        if te > data.time.max() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
        
        
        tslice=slice(ts, te)
        
        i0=np.abs(data.longitude.data-minlon).argmin() 
        i1=np.abs(data.longitude.data-maxlon).argmin() 
      
        j0=np.abs(data.latitude.data-minlat).argmin() 
        j1=np.abs(data.latitude.data-maxlat).argmin() 
        
        
        # expand the window a little bit        
        lon_0 = max(0, i0 - 2)
        lon_1 = min(data.longitude.size, i1 + 2)
        
        lat_0 = max(0, j0 - 2)
        lat_1 = min(data.latitude.size, j1 + 2)
                        
        if i0 > i1 :

            sh = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,data.longitude.size),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
            sh.longitude.values = sh.longitude.values -360.

            sh1 = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
            tot = xr.concat([sh,sh1],dim='longitude')
          
        else:            

            tot = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
#        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#            c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#            tot.longitude = tot.longitude + c*360.
      
            
        self.uvp = tot
      
      
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
        
        
        
    
class grib_pynio(meteo):
    
    def __init__(self,**kwargs):
    
    
        filenames = kwargs.get('mpaths', {})
        ft1 = kwargs.get('ft1', None)
        ft2 = kwargs.get('ft2', None)
        dft = kwargs.get('dft', None)
        engine = kwargs.get('engine', 'pynio')
        backend_kwargs = kwargs.get('backend_kwargs', {})
        
        kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['engine','backend_kwargs','mpaths','ft1','ft2','dft']}
        
        start_date = kwargs.get('start_date', None)
        try:
            start_date = pd.to_datetime(start_date)
        except:
            pass
    
        if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
        else:
          end_date = kwargs.get('end_date', None)
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass
        
        ts = pd.to_datetime(start_date)
        te = pd.to_datetime(end_date)     
        
                          
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
      
      #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo')
        sys.stdout.flush()
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, backend_kwargs=backend_kwargs, **kargs)    
        
        time_coord = data[data.data_vars.keys()[0]].dims[0]
        lon_coord = [x for (x,y) in data.dims.iteritems() if 'lon' in x][0]
        lat_coord = [x for (x,y) in data.dims.iteritems() if 'lat' in x][0]
                        
        u10 = [x for (x,y) in data.data_vars.iteritems() if ('10' in x) & ('u' in x.lower())][0]
        v10 = [x for (x,y) in data.data_vars.iteritems() if ('10' in x) & ('v' in x.lower())][0]
        msl = [x for (x,y) in data.data_vars.iteritems() if 'msl' in x.lower()][0]
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim=time_coord, backend_kwargs=backend_kwargs, **kargs)
        
        data = data.rename({msl:'msl',u10:'u10',v10:'v10',lon_coord:'longitude',lat_coord:'latitude',time_coord:'time'})
        
        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
        
                
        if not minlon : minlon = data.longitude.data.min()
        if not maxlon : maxlon = data.longitude.data.max()
        if not minlat : minlat = data.latitude.data.min()
        if not maxlat : maxlat = data.latitude.data.max()
        
        
        if minlon < data.longitude.data.min() : minlon = minlon + 360.
      
        if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.
      
        if minlon > data.longitude.data.max() : minlon = minlon - 360.
      
        if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
        
           
        if 'forecast' in time_coord:
            tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(data.time)
            data = data.assign_coords(time=tts)
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        
        
        if ts < data.time.min() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
          
        if te > data.time.max() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
        
        
        tslice=slice(ts, te)
        
        i0=np.abs(data.longitude.data-minlon).argmin() 
        i1=np.abs(data.longitude.data-maxlon).argmin() 
      
        j0=np.abs(data.latitude.data-minlat).argmin()
        j1=np.abs(data.latitude.data-maxlat).argmin() 
        
        # expand the window a little bit        
        lon_0 = max(0, i0 - 2)
        lon_1 = min(data.longitude.size, i1 + 2)
        
        lat_0 = max(0, j0 - 2)
        lat_1 = min(data.latitude.size, j1 + 2)
                
        if i0 > i1 :

            sh = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,data.longitude.size),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
            sh.longitude.values = sh.longitude.values -360.

            sh1 = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
            tot = xr.concat([sh,sh1],dim='longitude')
          
        else:            

            tot = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
#        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#            c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#            tot.longitude = tot.longitude + c*360.
      
            
        self.uvp = tot
      
      
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
        

def fix_my_data(ds):
    if 'A_PCP_GDS0_SFC_dif3h' not in ds.data_vars:
        mn = np.zeros(ds['PRMSL_GDS0_MSL'].shape)*np.nan
        ds['A_PCP_GDS0_SFC_dif3h'] = xr.DataArray(mn,coords=ds['PRMSL_GDS0_MSL'].coords, dims=ds['PRMSL_GDS0_MSL'].dims)
    return ds



class grib_pynio_AM(meteo):
    
    
    def __init__(self,**kwargs):
    
    
        filenames = kwargs.get('mpaths', {})
        ft1 = kwargs.get('ft1', None)
        ft2 = kwargs.get('ft2', None)
        dft = kwargs.get('dft', None)
        engine = kwargs.get('engine', 'pynio')
        backend_kwargs = kwargs.get('backend_kwargs', {})
        
        kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['engine','backend_kwargs','mpaths','ft1','ft2','dft']}
        
        start_date = kwargs.get('start_date', None)
        try:
            start_date = pd.to_datetime(start_date)
        except:
            pass
    
        if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
        else:
          end_date = kwargs.get('end_date', None)
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass
        
        ts = pd.to_datetime(start_date)
        te = pd.to_datetime(end_date)     
        
                          
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
      
      #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo')
        sys.stdout.flush()
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim='time', backend_kwargs=backend_kwargs, preprocess=fix_my_data,  **kargs)
        
        data = data.rename({'PRMSL_GDS0_MSL':'msl','U_GRD_GDS0_HTGL':'u10','V_GRD_GDS0_HTGL':'v10','g0_lon_1':'longitude','g0_lat_0':'latitude'})
        
        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
        
                
        if not minlon : minlon = data.longitude.data.min()
        if not maxlon : maxlon = data.longitude.data.max()
        if not minlat : minlat = data.latitude.data.min()
        if not maxlat : maxlat = data.latitude.data.max()
        
        
        if minlon < data.longitude.data.min() : minlon = minlon + 360.
      
        if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.
      
        if minlon > data.longitude.data.max() : minlon = minlon - 360.
      
        if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
        
        
        tt = [int(x.split('_')[-1]) for x in filenames]   
        tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(tt,unit='H')
        data = data.assign_coords(time=tts)
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        
        
        if ts < data.time.min() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
          
        if te > data.time.max() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
        
        
        tslice=slice(ts, te)
        
        i0=np.abs(data.longitude.data-minlon).argmin() 
        i1=np.abs(data.longitude.data-maxlon).argmin() 
      
        j0=np.abs(data.latitude.data-minlat).argmin()
        j1=np.abs(data.latitude.data-maxlat).argmin() 
        
        # expand the window a little bit        
        lon_0 = max(0, i0 - 2)
        lon_1 = min(data.longitude.size, i1 + 2)
        
        lat_0 = max(0, j0 - 2)
        lat_1 = min(data.latitude.size, j1 + 2)
                
        if i0 > i1 :

            sh = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,data.longitude.size),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
            sh.longitude.values = sh.longitude.values -360.

            sh1 = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
            tot = xr.concat([sh,sh1],dim='longitude')
          
        else:            

            tot = (
                data[['msl','u10', 'v10']]
                .isel(longitude=slice(lon_0,lon_1),latitude=slice(lat_0,lat_1))
                .sel(time=tslice)
                )
              
#        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#            c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#            tot.longitude = tot.longitude + c*360.
      
            
        self.uvp = tot
      
      
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
        


class grib_pynio_HNMS(meteo):
    
    
    def __init__(self,**kwargs):
    
    
        filenames = kwargs.get('mpaths', {})
        ft1 = kwargs.get('ft1', None)
        ft2 = kwargs.get('ft2', None)
        dft = kwargs.get('dft', None)
        engine = kwargs.get('engine', 'pynio')
        backend_kwargs = kwargs.get('backend_kwargs', {})
        
        kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['engine','backend_kwargs','mpaths','ft1','ft2','dft']}
        
        start_date = kwargs.get('start_date', None)
        try:
            start_date = pd.to_datetime(start_date)
        except:
            pass
    
        if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
        else:
          end_date = kwargs.get('end_date', None)
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass
        
        ts = pd.to_datetime(start_date)
        te = pd.to_datetime(end_date)     
        
                          
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
      
      #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo')
        sys.stdout.flush()
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim='time', backend_kwargs=backend_kwargs, **kargs)
        
        data = data.rename({'PS_MSL_GDS10_MSL':'msl','U_GDS10_HTGL':'u10','V_GDS10_HTGL':'v10','g10_lon_1':'longitude','g10_lat_0':'latitude'})
            
                
        if not minlon : minlon = data.longitude.data.min()
        if not maxlon : maxlon = data.longitude.data.max()
        if not minlat : minlat = data.latitude.data.min()
        if not maxlat : maxlat = data.latitude.data.max()
        
        
        if minlon < data.longitude.data.min() : minlon = minlon + 360.
      
        if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.
      
        if minlon > data.longitude.data.max() : minlon = minlon - 360.
      
        if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
        
        
        days = [int(x.split('lf')[1][:2]) for x in filenames]  
        hours = [int(x.split('lf')[1][2:4]) for x in filenames]
           
        tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(days,unit='D') + pd.to_timedelta(hours,unit='H')
        data = data.assign_coords(time=tts)
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        
        
        if ts < data.time.min() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
          
        if te > data.time.max() :
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('time frame not available\n')
            sys.stdout.write('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.stdout.flush()
            sys.exit(1)
        
        
        tslice=slice(ts, te)
        
        d1 = data.where(data.longitude>minlon,drop=True)
        d2 = d1.where(d1.longitude<maxlon,drop=True)
        d3 = d2.where(d2.latitude>minlat,drop=True)
        tot = d3.where(d3.latitude<maxlat,drop=True)

#        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#            c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#            tot.longitude = tot.longitude + c*360.
      
            
        self.uvp = tot[['msl','u10','v10','longitude','latitude']]      
      
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
        

    
class ecmwf_oper(meteo):   
        
    def __init__(self,**kwargs):
    
      filenames = kwargs.get('mpaths', {})
      ft1 = kwargs.get('ft1', None)
      ft2 = kwargs.get('ft2', None)
      dft = kwargs.get('dft', None)

      start_date = kwargs.get('start_date', None)
      try:
        start_date = pd.to_datetime(start_date)
      except:
        pass
    
      if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
      else:
          end_date = kwargs.get('end_date', None)
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass

                                              
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)   
      maxlat = kwargs.get('maxlat', None) 
          
      # read grib file and append to xarray

      pt=[]
      ut=[]
      vt=[]
      tt=[]


      sR = 0 # copmute the time step of the data
      dsR = 0

      
      for filename in filenames:
        
        #--------------------------------------------------------------------- 
#        sys.stdout.flush()
#        sys.stdout.write('\n')
#        sys.stdout.write('extracting meteo from {}\n'.format(filename))
#        sys.stdout.flush()
        #---------------------------------------------------------------------        

     	try: 
       		f = open(filename)
      	except:
        	print('no file {}'.format(filename))
         	sys.exit(1)
        while True:
          try:
            gid = grib_new_from_file(f)#,headers_only = True)
            if gid is None: 
                sys.stdout.write('end of file {}\n'.format(filename))
                break
            
            date=grib_get(gid, 'date')
            dataTime=grib_get(gid, 'dataTime')
            stepRange=grib_get(gid, 'stepRange')
            dsR_=int(stepRange)-sR # current time step
            dsR = max(dsR_, dsR) #save the max 
            sR=int(stepRange) #update for future comparison
            timestamp = pd.to_datetime(str(date)) + pd.to_timedelta('{}H'.format(dataTime/100.))
            tstamp = timestamp+pd.to_timedelta('{}H'.format(stepRange))
            try:
                if (ft1 <= int(stepRange) <= ft2) & (tstamp <= self.end_date): #+ pd.to_timedelta('{}H'.format(dsR))):
                
                    name,varin,ilon,ilat=getd(gid)    
            
                else:
                    grib_release(gid)
                    if tstamp > self.end_date : break #+ pd.to_timedelta('{}H'.format(dsR)): break
                    continue
            except:
                name,varin,ilon,ilat=getd(gid)
                
            lon=ilon[0,:]
            lat=ilat[:,0]
           
            if not minlon : minlon = lon.min()
            if not maxlon : maxlon = lon.max()
            if not minlat : minlat = lat.min()
            if not maxlat : maxlat = lat.max()
            
        # verbose 
        #--------------------------------------------------------------------- 
#            sys.stdout.flush()
#            sys.stdout.write('\n')
#            sys.stdout.write('retrieving {} at {}\n'.format(name, tstamp))
#            sys.stdout.flush()
        #---------------------------------------------------------------------        
        
        
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
                     tt.append(tstamp)
            elif name == '10u':
                     ut.append(data)
            elif name == '10v':
                     vt.append(data)


    # END OF FOR

          except Exception as e:
            print(e)
            print('ERROR in meteo file {}'.format(date))

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
    
      self.uvp = met
      self.ft1 = ft1
      self.ft2 = ft2 
      self.dft = dft   
      
      self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

      self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      #--------------------------------------------------------------------- 
#      sys.stdout.flush()
#      sys.stdout.write('\n')
#      sys.stdout.write('meteo done\n')
#      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
    
        
            
      

class gfs_erdap(meteo):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
      ts = kwargs.get('start_date', None)
      te = kwargs.get('end_date', None)
      
      ts = pd.to_datetime(ts)
      te = pd.to_datetime(te)     
         
      url = kwargs.get('meteo_url', 'https://bluehub.jrc.ec.europa.eu/erddap/griddap/NCEP_Global_Best')
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting meteo from {}.html\n'.format(url))
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
              .isel(longitude=slice(i0,data.longitude.size),latitude=slice(j0,j1+1))
              .sel(time=tslice)
              )
          sh.longitude.values = sh.longitude.values -360.

          sh1 = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(longitude=slice(0,i1+1),latitude=slice(j0,j1+1))
              .sel(time=tslice)
              )
              
          tot = xr.concat([sh,sh1],dim='longitude')
          
      else:            

          tot = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(longitude=slice(i0,i1+1),latitude=slice(j0,j1+1))
              .sel(time=tslice)
              )
              
      if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
          c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
          tot['longitude'] = tot['longitude'] + c*360.
      
      xx,yy=np.meshgrid(tot.longitude,tot.latitude)
      
      other = xr.Dataset({'lons': (['x', 'y'], xx), 
                          'lats': (['x', 'y'], yy) })
      
      self.uvp = xr.merge([tot, other])
      
      self.uvp.attrs =  tot.attrs
      
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
    
        filename = kwargs.get('mpaths', None)
        
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo from {}\n'.format(filename))
        sys.stdout.flush()
        #---------------------------------------------------------------------      
        
        
        self.uvp = xr.open_dataset(filename)
         
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
        
        
        
class hnms_oper(meteo):
    
    def __init__(self,**kwargs):
    
        filenames = kwargs.get('mpaths', {})
        
        ncores = kwargs.get('ncores', 1)
    
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)
        maxlat = kwargs.get('maxlat', None) 
        
        try:
            minlon = minlon - .1
            maxlon = maxlon + .1
            minlat = minlat - .1
            maxlat = maxlat + .1
        except:
            pass
  
        # read the first file to get variables
        f = open(filenames[0])
        gid = grib_new_from_file(f)
        
        
        Ni=grib_get(gid,'Ni')
        Nj=grib_get(gid,'Nj')
                
        elat=grib_get_array(gid,'latitudes')
        elon=grib_get_array(gid,'longitudes')
        
        elon = elon.reshape(Nj,Ni)
        elat = elat.reshape(Nj,Ni)
        
        orig = geometry.SwathDefinition(lons=elon,lats=elat) # original grid
        
        gridType = grib_get(gid,'gridType')
        
        if not minlon : minlon = elon.min()
        if not maxlon : maxlon = elon.max()
        if not minlat : minlat = elat.min()
        if not maxlat : maxlat = elat.max()
        
        
        
        #Set lat/lon window for interpolation
        prj = pyproj.Proj('+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
        [[a0,a1],[a2,a3]] = prj([minlon, minlat], [maxlon, maxlat])
  
        area_id = 'HNMS'
        description = 'HNMS COSMO'
        proj_id = 'eqc'
        projection = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
        x_size = Ni
        y_size = Nj
        area_extent = (a0, a1, a2, a3)
        target_def = utils.get_area_def(area_id, description, proj_id, projection,
                         x_size, y_size, area_extent)

        lons,lats = geometry.AreaDefinition.get_lonlats(target_def) # out grid
        #compute bilinear interoplation parameters
        t_params, s_params, input_idxs, idx_ref = \
                              bilinear.get_bil_info(orig, target_def, radius=50e3, nprocs=1)
  
        grib_release(gid)
        f.close()                     
  
        # read grib file and append to xarray
        pt=[]
        ut=[]
        vt=[]
        tt=[]

        for filename in filenames:

          #--------------------------------------------------------------------- 
          #  sys.stdout.flush()
          #  sys.stdout.write('\n')
          #  sys.stdout.write('extracting meteo from {}\n'.format(filename))
          #  sys.stdout.flush()
          #---------------------------------------------------------------------      


            try: 
                f = open(filename)
            except:
                print('no file {}'.format(filename))
                sys.exit(1)



            try:    
                for it in range(3):
        
                    gid = grib_new_from_file(f)
    
                    name=grib_get(gid, 'shortName')
                    mv=grib_get(gid,'missingValue')

                    date=grib_get(gid, 'date')
                    dataTime=grib_get(gid, 'dataTime')
                    stepRange=grib_get(gid, 'stepRange')

                    values=grib_get_values(gid)    

                    q = values.reshape(Nj,Ni)
                       
                    data = bilinear.get_sample_from_bil_info(q.ravel(), t_params, s_params,
                                                             input_idxs, idx_ref,
                                                             output_shape=target_def.shape)
    
#                    data = bilinear.resample_bilinear(q, orig, target_def,
#                                           radius=50e3, neighbours=32,
#                                           nprocs=ncores, fill_value=0,
#                                           reduce_data=False, segments=None,
#                                           epsilon=0)
    
                    timestamp =  datetime.datetime.strptime(str(date),'%Y%m%d')+datetime.timedelta(hours=dataTime/100) 


                    if name == 'msl' : 
                         pt.append(data)
                         tt.append(timestamp+datetime.timedelta(hours=int(stepRange)))
                    elif name == '10u':
                         ut.append(data)
                    elif name == '10v':
                         vt.append(data)


            # END OF FOR

            except Exception as e:
                print(e)
                print('ERROR in meteo file {}'.format(date))

            f.close()

        met = xr.Dataset({'msl': (['time', 'latitude', 'longitude'],  np.array(pt)), 
                              'u10': (['time', 'latitude', 'longitude'], np.array(ut)),   
                              'v10': (['time', 'latitude', 'longitude'], np.array(vt)),   
                              'lons': (['x', 'y'], lons),   
                              'lats': (['x', 'y'], lats)},   
                              coords={'longitude': ('longitude', lons[0,:]),   
                                      'latitude': ('latitude', lats[:,0]),   
                                      'time': tt })
        
        #mask non values
        met['msl'] = met.msl.where(met.msl>0)
        met['u10'] = met.u10.where(met.msl>0)
        met['v10'] = met.v10.where(met.msl>0)
                                      
        self.uvp = met
    
        self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

        self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      #--------------------------------------------------------------------- 
      #  sys.stdout.flush()
      #  sys.stdout.write('\n')
      #  sys.stdout.write('meteo done\n')
      #  sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
    

class am_oper(meteo):
    
    def __init__(self,**kwargs):
    
        filenames = kwargs.get('mpaths', {})
        
        ncores = kwargs.get('ncores', 1)
              
        # read the first file to get variables
        f = open(filenames[0])
        gid = grib_new_from_file(f)
        
        
        lonfgp=grib_get(gid,'longitudeOfFirstGridPointInDegrees')
        latfgp=grib_get(gid,'latitudeOfFirstGridPointInDegrees')
        lonlgp=grib_get(gid,'longitudeOfLastGridPointInDegrees')
        latlgp=grib_get(gid,'latitudeOfLastGridPointInDegrees')

        
        Ni=grib_get(gid,'Ni')
        Nj=grib_get(gid,'Nj')
                
        elat=grib_get_array(gid,'latitudes')
        elon=grib_get_array(gid,'longitudes')
        
        if latfgp > latlgp :
            elat = elat[::-1]
        
        lon = elon.reshape(Nj,Ni)
        lat = elat.reshape(Nj,Ni)
        
        gridType = grib_get(gid,'gridType')
        
        minlon = kwargs.get('minlon', elon.min())
        maxlon = kwargs.get('maxlon', elon.max())
        minlat = kwargs.get('minlat', elat.min())
        maxlat = kwargs.get('maxlat', elat.max()) 


        grib_release(gid)
        f.close()                     
  

  
        if (minlon < elon.min()) or (maxlon > elon.max()) or (minlat < elat.min()) or (maxlat > elat.max()): 
            print(minlon,maxlon,minlat,maxlat)
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('Meteo Problem!\n')
            sys.stdout.write('Lon must be within {} and {}\n'.format(elon.min(),elon.max()))
            sys.stdout.write('Lat must be within {} and {}\n'.format(elat.min(),elat.max()))
            sys.stdout.flush()
            sys.exit(1)
            
        # get range of values
        i1=np.abs(lon[0,:]-minlon).argmin()-2
        i2=np.abs(lon[0,:]-maxlon).argmin()+2
        j1=np.abs(lat[:,0]-minlat).argmin()-2
        j2=np.abs(lat[:,0]-maxlat).argmin()+2
    
        if i1 < 0 : i1 = 0 # fix limits
        if i2 > lon.shape[1] : i2 = lon.shape[1]   
        if j1 < 0 : j1 = 0
        if j2 > lat.shape[0]: j2 = lat.shape[0]
        
        # read grib file and append to xarray
        pt=[]
        ut=[]
        vt=[]
        tt=[]

        for filename in filenames:

        #--------------------------------------------------------------------- 
        #    sys.stdout.flush()
        #    sys.stdout.write('\n')
        #    sys.stdout.write('extracting meteo from {}\n'.format(filename))
        #    sys.stdout.flush()
        #---------------------------------------------------------------------      

            try: 
                f = open(filename)
            except:
                print('no file {}'.format(filename))
                sys.exit(1)



            try:    
                for it in range(3):
        
                    gid = grib_new_from_file(f)
    
                    name=grib_get(gid, 'shortName')
                    mv=grib_get(gid,'missingValue')

                    date=grib_get(gid, 'date')
                    dataTime=grib_get(gid, 'dataTime')
                    stepRange=grib_get(gid, 'stepRange')

                    values=grib_get_values(gid)    
                    
                    q = values.reshape(Nj,Ni)
                    
                    if latfgp > latlgp :
                        q=np.flipud(q)

                           
                    timestamp =  datetime.datetime.strptime(str(date),'%Y%m%d')+datetime.timedelta(hours=dataTime/100) 


                    lons, lats = lon[j1:j2,i1:i2],lat[j1:j2,i1:i2]
                    data = deepcopy(q[j1:j2,i1:i2])


                    if name == 'pmsl' : 
                         pt.append(data)
                         tt.append(timestamp+datetime.timedelta(hours=int(stepRange)))
                    elif name == '10u':
                         ut.append(data)
                    elif name == '10v':
                         vt.append(data)

            # END OF FOR

            except Exception as e:
                print(e)
                print('ERROR in meteo file {}'.format(date))

            f.close()

            
        met = xr.Dataset({'msl': (['time', 'latitude', 'longitude'],  np.array(pt)), 
                              'u10': (['time', 'latitude', 'longitude'], np.array(ut)),   
                              'v10': (['time', 'latitude', 'longitude'], np.array(vt)),   
                              'lons': (['x', 'y'], lons),   
                              'lats': (['x', 'y'], lats)},   
                              coords={'longitude': ('longitude', lons[0,:]),   
                                      'latitude': ('latitude', lats[:,0]),   
                                      'time': tt })
                                              
        self.uvp = met
    
        self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

        self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      #--------------------------------------------------------------------- 
      #  sys.stdout.flush()
      #  sys.stdout.write('\n')
      #  sys.stdout.write('meteo done\n')
      #  sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
 
 

class gfs_oper(meteo):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
      ts = kwargs.get('start_date', None)
      te = kwargs.get('end_date', None)
      
      ts = pd.to_datetime(ts)
      te = pd.to_datetime(te)     
      
      
      url0='http://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs{}/gfs_0p25_1hr_{}z'.format(ts.strftime('%Y%m%d'),ts.strftime('%H'))
      
      url = kwargs.get('meteo_url', url0)

      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting meteo from {}\n'.format(url))
      sys.stdout.flush()
      #---------------------------------------------------------------------      

      try:
         data = xr.open_dataset(url)    
      except:
         #--------------------------------------------------------------------- 
         sys.stdout.flush()
         sys.stdout.write('\n')
         sys.stdout.write('Please provide date data within the last 10 days\n')
         sys.stdout.flush()
         #---------------------------------------------------------------------      
         sys.exit(1)
          
      
      if minlon < data.lon.minimum : minlon = minlon + 360.
      
      if maxlon < data.lon.minimum : maxlon = maxlon + 360.
      
      if minlon > data.lon.maximum : minlon = minlon - 360.
      
      if maxlon > data.lon.maximum : maxlon = maxlon - 360.
                  
      tslice=slice(ts, te)
    
      i0=np.abs(data.lon.data-minlon).argmin()
      i1=np.abs(data.lon.data-maxlon).argmin()

      
      j0=np.abs(data.lat.data-minlat).argmin()
      j1=np.abs(data.lat.data-maxlat).argmin()

      if i0 > i1 :

          sh = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(lon=slice(i0,data.lon.size),lat=slice(j0,j1+1))
              .sel(time=tslice)
              )
          sh.lon.values = sh.lon.values -360.

          sh1 = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(lon=slice(0,i1+1),lat=slice(j0,j1+1))
              .sel(time=tslice)
              )
              
          tot = xr.concat([sh,sh1],dim='lon')
          
      else:            

          tot = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(lon=slice(i0,i1+1),lat=slice(j0,j1+1))
              .sel(time=tslice)
              )
              
      if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
          c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
          tot['lon'] = tot['lon'] + c*360.
      
      xx,yy=np.meshgrid(tot.lon,tot.lat)
      
      other = xr.Dataset({'lons': (['x', 'y'], xx), 
                          'lats': (['x', 'y'], yy) })
      
      self.uvp = xr.merge([tot, other])
      
      self.uvp.rename({'prmslmsl':'msl','ugrd10m':'u10','vgrd10m':'v10'}, inplace=True)
      
      self.hview = hv.Dataset(self.uvp,kdims=['time','lon','lat'],vdims=['msl','u10','v10'])

      self.gview = gv.Dataset(self.uvp,kdims=['time','lon','lat'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      
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
        
    

class erai(meteo):
    
    def __init__(self,**kwargs):
        
      filenames = kwargs.get('mpaths', {})
      ft1 = kwargs.get('ft1', None)
      ft2 = kwargs.get('ft2', None)
      dft = kwargs.get('dft', None)
                  
      start_date = kwargs.get('start_date', None)
      if not start_date : start_date = '1970-01-01'
      try:
            start_date = pd.to_datetime(start_date)
      except:
            pass
    
      if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          end_date = start_date + pd.to_timedelta(time_frame)
      else:
          end_date = kwargs.get('end_date', None)
          if not end_date : end_date = '2070-01-01'
          try:
              end_date = pd.to_datetime(end_date)
              time_frame = end_date - start_date
          except:
              pass
      
      
                                      
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)   
      maxlat = kwargs.get('maxlat', None) 
              
       # read grib file and append to xarray

      pt=[]
      ut=[]
      vt=[]
      tt=[]

          
      for filename in filenames:
        
          #--------------------------------------------------------------------- 
#          sys.stdout.flush()
#          sys.stdout.write('\n')
#          sys.stdout.write('extracting meteo from {}\n'.format(filename))
#          sys.stdout.flush()
          #---------------------------------------------------------------------        

          try: 
              f = open(filename)
          except:
              print('no file {}'.format(filename))
              sys.exit(1)
          while True:
              try:
                  gid = grib_new_from_file(f)#,headers_only = True)
                  if gid is None: 
                      sys.stdout.write('end of file {}\n'.format(filename))
                      break
            
                  date=grib_get(gid, 'date')
                  dataTime=grib_get(gid, 'dataTime')
                  stepRange=grib_get(gid, 'stepRange')
                  timestamp = pd.to_datetime(str(date)) + pd.to_timedelta('{}H'.format(dataTime/100.))
                  tstamp = timestamp+pd.to_timedelta('{}H'.format(stepRange))
                  if (tstamp >= pd.to_datetime(start_date)) & (tstamp <= pd.to_datetime(end_date)):
                
                      name,varin,ilon,ilat=getd(gid)    

                  else:
                      grib_release(gid)
                      if tstamp > pd.to_datetime(end_date) : break
                      continue
                
                  lon=ilon[0,:]
                  lat=ilat[:,0]
           
              # verbose 
              #--------------------------------------------------------------------- 
#                  sys.stdout.flush()
#                  sys.stdout.write('\n')
#                  sys.stdout.write('retrieving {} at {}\n'.format(name, tstamp))
#                  sys.stdout.flush()
              #--------------------------------------------------------------------- 
                     
                  if not minlon : minlon = lon[0]
                  if not maxlon : maxlon = lon[-1]
                  if not minlat : minlat = lat[0]
                  if not maxlat : maxlat = lat[-1]
                     
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
                           tt.append(tstamp)
                  elif name == '10u':
                           ut.append(data)
                  elif name == '10v':
                           vt.append(data)


          # END OF FOR

              except Exception as e:
                  print(e)
                  print('ERROR in meteo file {}'.format(date))

          f.close()

      met = xr.Dataset({'msl': (['time', 'latitude', 'longitude'],  np.array(pt)), 
                                'u10': (['time', 'latitude', 'longitude'], np.array(ut)),   
                                'v10': (['time', 'latitude', 'longitude'], np.array(vt)),   
                                'lons': (['x', 'y'], lons),   
                                'lats': (['x', 'y'], lats)},   
                                coords={'longitude': ('longitude', lons[0,:]),   
                                        'latitude': ('latitude', lats[:,0]),   
                                        'time': tt })   

      self.uvp = met
      self.ft1 = ft1
      self.ft2 = ft2 
      self.dft = dft   
      
      self.hview = hv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'])

      self.gview = gv.Dataset(self.uvp,kdims=['time','longitude','latitude'],vdims=['msl','u10','v10'],crs=crs.PlateCarree())
      
      
      #--------------------------------------------------------------------- 
#      sys.stdout.flush()
#      sys.stdout.write('\n')
#      sys.stdout.write('meteo done\n')
#      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
    
        
      
    