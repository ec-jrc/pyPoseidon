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
import os
import time
import xarray as xr
import pandas as pd
import importlib
from pyPoseidon.utils.get_value import get_value
from dateutil import parser
import xesmf as xe
import logging

#logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(name)s:%(message)s')

file_handler = logging.FileHandler('meteo.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

sformatter = logging.Formatter('%(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(sformatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)

def fix_my_data(ds):
    if 'A_PCP_GDS0_SFC_dif3h' not in ds.data_vars:
        mn = np.zeros(ds['PRMSL_GDS0_MSL'].shape)*np.nan
        ds['A_PCP_GDS0_SFC_dif3h'] = xr.DataArray(mn,coords=ds['PRMSL_GDS0_MSL'].coords, dims=ds['PRMSL_GDS0_MSL'].dims)
    return ds



def regrid(ds):
    

    de = ds.rename({'longitude': 'lon', 'latitude': 'lat'})

    minlon = de.lon.values.min() # get lat/lon window
    maxlon = de.lon.values.max()
    minlat = de.lat.values.min()
    maxlat = de.lat.values.max()

    Nj, Ni = de.lon.shape # get shape of original grid
    
    y = np.linspace(minlat, maxlat, Nj) # create a similar grid as the original
    x = np.linspace(minlon, maxlon, Ni)
    
    dlon = np.diff(x)[0] #resolution
    dlat = np.diff(y)[0]
    
    ds_out = xe.util.grid_2d(minlon,maxlon,dlon,minlat,maxlat,dlat) # establish out grid 
    
    regridder = xe.Regridder(de, ds_out, 'bilinear') 
    
    varn = ['msl','u10','v10']

    dats = []
    for var in varn:
         
            q = regridder(de[var])
          
            dats.append(q.where(q.values != 0))
    
    data = xr.merge(dats)
    
    data = data.assign_coords(x=data.lon.values[0,:],y=data.lat.values[:,0])# assign 1-D coords
    
    data = data.rename({'x':'longitude','y':'latitude'}).drop(['lon','lat']) # rename and drop 2-D coords
    
    return data



class meteo:
    impl=None
    def __init__(self,**kwargs):
        msource = kwargs.get('meteo', None)
        combine = kwargs.get('combine', False)
        if combine :
            self.impl = combine_meteo(**kwargs)
        else:
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
        engine = kwargs.get('engine', None) 
        service = kwargs.get('service', None)   

        if engine == 'cfgrib' :  
                         
            dats = []
            for filename in filenames[:-1]:
            
                kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
                kargs['mpaths'] = [filename]
                
                ds = grib_cfgrib(**kargs)
                dats.append(ds.uvp)
        
            kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
            kargs['ft2'] += 1
            kargs['mpaths'] = [filenames[-1]]
            ds = grib_cfgrib(**kargs)
            dats.append(ds.uvp)
                
                
        elif engine == 'pynio' :
                   
            
          if service == 'HNMS' :
            
            dats = []
            for filename in filenames[:-1]:
            
                kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
                kargs['mpaths'] = [filename]
                
                ds = grib_pynio_HNMS(**kargs)
                dats.append(ds.uvp)
        
            kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
            kargs['ft2'] += 1
            kargs['mpaths'] = [filenames[-1]]
            
            ds = grib_pynio_HNMS(**kargs)
            dats.append(ds.uvp)    
        
        
          elif service == 'AM' :

            dats = []
            for filename in filenames[:-1]:
            
                kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
                kargs['mpaths'] = [filename]
                
                ds = grib_pynio_AM(**kargs)
                dats.append(ds.uvp)
        
            kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
            kargs['ft2'] += 1
            kargs['mpaths'] = [filenames[-1]]
            
            ds = grib_pynio_AM(**kargs)
            dats.append(ds.uvp)
            
          else:    
            
            dats = []
            for filename in filenames[:-1]:
            
                kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
                kargs['mpaths'] = [filename]
                
                ds = grib_pynio(**kargs)
                dats.append(ds.uvp)
        
            kargs = { key: kwargs[key] for key in kwargs.keys() if key not in ['start_date','end_date']}
            kargs['ft2'] += 1
            kargs['mpaths'] = [filenames[-1]]
            
            ds = grib_pynio(**kargs)
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
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
        xr_kwargs = kwargs.get('xr_kwargs', {'concat_dim':'step'})
    
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
      
      #--------------------------------------------------------------------- 
        logger.info('extracting meteo')
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, backend_kwargs=backend_kwargs, **xr_kwargs)    
        
        data = data.squeeze(drop=True)
#        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing> not efficient for output      
        
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
                            
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        if ft1 :
            ts = data.time.data[ft1]
        
        if ft2 :
            te = data.time.data[ft2]
        
        
        if ts < data.time.data.min() :
            logger.warning('coverage between {} and {} \n'.format(pd.to_datetime(data.valid_time.min().values).strftime('%Y.%m.%d %H:%M:%S'),pd.to_datetime(data.valid_time.max().values).strftime('%Y.%m.%d %H:%M:%S')))
            logger.error('time frame not available')
            sys.exit(1)
          
        if te > data.time.data.max() :
            logger.warning('coverage between {} and {} \n'.format(pd.to_datetime(data.valid_time.min().values).strftime('%Y.%m.%d %H:%M:%S'),pd.to_datetime(data.valid_time.max().values).strftime('%Y.%m.%d %H:%M:%S')))
            logger.error('time frame not available\n')
            sys.exit(1)
        
        
        if len(data.longitude.shape) == 2:
            d1 = data.where(data.longitude>minlon,drop=True)
            d2 = d1.where(d1.longitude<maxlon,drop=True)
            d3 = d2.where(d2.latitude>minlat,drop=True)
            d4 = d3.where(d3.latitude<maxlat,drop=True)        
            
            data = regrid(d4)
                    
        
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
                   
        # descenting lats
        if j0 > j1 : 
           j0, j1 = j1, j0
           lat_0 = max(0, j0 - 1)
           lat_1 = min(data.latitude.size, j1 + 3)
           
         
        
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
        logger.info('meteo done\n')
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
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
        xr_kwargs = kwargs.get('xr_kwargs', {})
        if 'preprocess' in xr_kwargs.keys():
            xr_kwargs['preprocess'] = fix_my_data
        print(xr_kwargs)
        
        
        
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
      
      #--------------------------------------------------------------------- 
        logger.info('extracting meteo')
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, backend_kwargs=backend_kwargs, **xr_kwargs)    
        
        key = [key for key in data.data_vars.keys()][0]
        time_coord = data[key].dims[0]
        lon_coord = [x for (x,y) in data.dims.items() if 'lon' in x][0]
        lat_coord = [x for (x,y) in data.dims.items() if 'lat' in x][0]
                        
        u10 = [x for (x,y) in data.data_vars.items() if ('10' in x) | ('u' in x.lower())][0]
        v10 = [x for (x,y) in data.data_vars.items() if ('10' in x) | ('v' in x.lower())][0]
        msl = [x for (x,y) in data.data_vars.items() if 'msl' in x.lower()][0]
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim=time_coord, backend_kwargs=backend_kwargs, **xr_kwargs)
        
        data = data.rename({msl:'msl',u10:'u10',v10:'v10',lon_coord:'longitude',lat_coord:'latitude',time_coord:'time'})
        
#        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
        
                
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
        
        if ft1 :
            ts = data.time.data[ft1]
        
        if ft2 :
            te = data.time.data[ft2]
        
        
        if ts < data.time.data.min() :
            logger.warning('coverage between {} and {} \n'.format(data.time.data.min(),data.time.data.max()))
            logger.error('time frame not available\n')
            sys.exit(1)
          
        if te > data.time.data.max() :
            logger.error('time frame not available\n')
            logger.warning('coverage between {} and {} \n'.format(data.time.data.min(),data.time.data.max()))
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
        
        # descenting lats
        if j0 > j1 : 
           j0, j1 = j1, j0
           lat_0 = max(0, j0 - 1)
           lat_1 = min(data.latitude.size, j1 + 3)
           
        
                
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
        logger.info('meteo done\n')
        #--------------------------------------------------------------------- 
        
        
        
    def output(self,solver=None,**kwargs):
         
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
        



class grib_pynio_AM(meteo):
    
    
    def __init__(self,**kwargs):
    
    
        filenames = kwargs.get('mpaths', {})
        ft1 = kwargs.get('ft1', None)
        ft2 = kwargs.get('ft2', None)
        dft = kwargs.get('dft', None)
        engine = kwargs.get('engine', 'pynio')
        backend_kwargs = kwargs.get('backend_kwargs', {})
        
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
        xr_kwargs = kwargs.get('xr_kwargs', {})
              
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
        
                          
      #--------------------------------------------------------------------- 
        logger.info('extracting meteo')
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim='time', backend_kwargs=backend_kwargs, preprocess=fix_my_data,  **xr_kwargs)
        
        data = data.rename({'PRMSL_GDS0_MSL':'msl','U_GRD_GDS0_HTGL':'u10','V_GRD_GDS0_HTGL':'v10','g0_lon_1':'longitude','g0_lat_0':'latitude'})
        
#        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
        
                
        if not minlon : minlon = data.longitude.data.min()
        if not maxlon : maxlon = data.longitude.data.max()
        if not minlat : minlat = data.latitude.data.min()
        if not maxlat : maxlat = data.latitude.data.max()
        
        
#        if minlon < data.longitude.data.min() : minlon = minlon + 360.
#      
#        if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.
#      
#        if minlon > data.longitude.data.max() : minlon = minlon - 360.
#      
#        if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
        
        print(data.msl.attrs)
        tt = [int(x.split('_')[-1]) for x in filenames]   
        tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(tt,unit='H')
        data = data.assign_coords(time=tts)
            
        if not ts : ts = data.time.data.min()
        if not te : te = data.time.data.max()
        
        
        
        if ts < data.time.data.min() :
            logger.warning('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
            logger.error('time frame not available\n')
            sys.exit(1)
          
        if te > data.time.data.max() :
            logger.error('time frame not available\n')
            logger.warning('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
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
        
        # descenting lats
        if j0 > j1 : 
           j0, j1 = j1, j0
           lat_0 = max(0, j0 - 1)
           lat_1 = min(data.latitude.size, j1 + 3)
                           
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
        logger.info('meteo done\n')
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
        
        minlon = kwargs.get('minlon', None)
        maxlon = kwargs.get('maxlon', None)
        minlat = kwargs.get('minlat', None)   
        maxlat = kwargs.get('maxlat', None) 
        xr_kwargs = kwargs.get('xr_kwargs', {})
                
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
        
                               
      #--------------------------------------------------------------------- 
        logger.info('extracting meteo')
      #---------------------------------------------------------------------      
        
        data = xr.open_mfdataset(filenames, engine=engine, concat_dim='time', backend_kwargs=backend_kwargs, **xr_kwargs)
        
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
        
        
        
        if ts < data.time.data.min() :
            logger.error('time frame not available\n')
            logger.warning('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.exit(1)
          
        if te > data.time.data.max() :
            logger.error('time frame not available\n')
            logger.warning('coverage between {} and {} \n'.format(data.valid_time.min(),data.time.max()))
            sys.exit(1)
        
        
        d1 = data.where(data.longitude>minlon,drop=True)
        d2 = d1.where(d1.longitude<maxlon,drop=True)
        d3 = d2.where(d2.latitude>minlat,drop=True)
        tot = d3.where(d3.latitude<maxlat,drop=True)        

#        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#            c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#            tot.longitude = tot.longitude + c*360.
      
        elon = tot.longitude.values
        elat = tot.latitude.values
        
        Ni, Nj = tot.longitude.shape
        
        orig = geometry.SwathDefinition(lons=elon,lats=elat) # original grid
          
       #Set lat/lon window for interpolation
        prj = pyproj.Proj('+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs')
        [[a0,a1],[a2,a3]] = prj([minlon, minlat], [maxlon, maxlat])

        orig = geometry.SwathDefinition(lons=elon,lats=elat) # original points


        area_id = 'HNMS'
        description = 'HNMS COSMO'
        proj_id = 'eqc'
        projection = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
        x_size = Ni
        y_size = Nj
        area_extent = (a0, a1, a2, a3)
        target_def = get_area_def(area_id, description, proj_id, projection,
                                   x_size, y_size, area_extent)

        lons,lats = geometry.AreaDefinition.get_lonlats(target_def)

        t_params, s_params, input_idxs, idx_ref = \
                              bilinear.get_bil_info(orig, target_def, radius=50e3, nprocs=1)

        varn = ['msl','u10','v10']
  
        for var in varn:
            
                    q = tot[var].values
    
                    newdata = bilinear.get_sample_from_bil_info(q.ravel(), t_params, s_params,
                                                 input_idxs, idx_ref,
                                                 output_shape=target_def.shape)

                    me = xr.DataArray(newdata, coords={'latitude': lats[:,0], 'longitude': lons[0,:]},
                     dims=['latitude', 'longitude'])
            
                    tot[var] = me   
            
        
        self.uvp = tot.drop('g10_rot_2')      
      
       #--------------------------------------------------------------------- 
        logger.info('meteo done\n')
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
         
      url = kwargs.get('meteo_url', 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best')
      
      #--------------------------------------------------------------------- 
      logger.info('extracting meteo from {}.html\n'.format(url))
      #---------------------------------------------------------------------      

      data = xr.open_dataset(url)    
      
      if minlon < data.geospatial_lon_min : minlon = minlon + 360.
      
      if maxlon < data.geospatial_lon_min : maxlon = maxlon + 360.
      
      if minlon > data.geospatial_lon_max : minlon = minlon - 360.
      
      if maxlon > data.geospatial_lon_max : maxlon = maxlon - 360.
            
      if ts < parser.parse(data.attrs['time_coverage_start']).replace(tzinfo=None) :
          logger.error('time frame not available\n')
          logger.warning('coverage between {} and {} \n'.format(data.attrs['time_coverage_start'],data.attrs['time_coverage_end']))
          sys.exit(1)
          
      if te > parser.parse(data.attrs['time_coverage_end']).replace(tzinfo=None) :
          logger.error('time frame not available\n')
          logger.warning('coverage between {} and {} \n'.format(data.attrs['time_coverage_start'],data.attrs['time_coverage_end']))
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
      
      self.uvp = tot
      
      self.uvp.attrs =  tot.attrs
      
      self.uvp = self.uvp.rename({'prmslmsl':'msl','ugrd10m':'u10','vgrd10m':'v10'})
            
      
      
      #--------------------------------------------------------------------- 
      logger.info('meteo done\n')
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
        logger.info('extracting meteo from {}\n'.format(filename))
        #---------------------------------------------------------------------      
        
        
        self.uvp = xr.open_dataset(filename)
         
      #--------------------------------------------------------------------- 
        logger.info('meteo done\n')
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
      
      
      url0='https://nomads.ncep.noaa.gov:9090/dods/gfs_0p25_1hr/gfs{}/gfs_0p25_1hr_{}z'.format(ts.strftime('%Y%m%d'),ts.strftime('%H'))
      
      url = kwargs.get('meteo_url', url0)

      #--------------------------------------------------------------------- 
      logger.info('extracting meteo from {}\n'.format(url))
      #---------------------------------------------------------------------      

      try:
         data = xr.open_dataset(url)    
      except:
         #--------------------------------------------------------------------- 
         logger.error('Please provide date data within the last 10 days\n')
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
            
      self.uvp = tot
      
      self.uvp = self.uvp.rename({'prmslmsl':'msl','ugrd10m':'u10','vgrd10m':'v10','lon':'longitude','lat':'latitude'})
      
      
      #--------------------------------------------------------------------- 
      logger.info('meteo done\n')
      #--------------------------------------------------------------------- 
      
      
    def output(self,solver=None,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
                
        model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
        s = getattr(model,solver) # get solver class
                
        s.to_force(self.uvp,vars=['msl','u10','v10'],rpath=path)
        
    

