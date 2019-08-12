# -*- coding: utf-8 -*- 
"""
Meteo module. Pre-processing the weather forcing component.

"""

# Copyright 2018 European Union
# This file is part of pyPoseidon.
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
import xesmf as xe
import logging

#logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(levelname)-8s %(asctime)s:%(name)s:%(message)s')

file_handler = logging.FileHandler('meteo.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

sformatter = logging.Formatter('%(levelname)-8s %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(sformatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)



def erddap():
    
    links = [
            'http://oos.soest.hawaii.edu/erddap/griddap/NCEP_Global_Best',
            'https://cidportal.jrc.ec.europa.eu/services/erddap/griddap/NCEP_Global_Best',
            'https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best',
             ]
    return links


def reduced_gg(dc):
    
    logger.info('regriding meteo')
    
    dc = dc.sortby('latitude', ascending=True)
    
    
    minlon = dc.longitude.values.min()
    maxlon = dc.longitude.values.max()
    minlat = dc.latitude.values.min()
    maxlat = dc.latitude.values.max()
    
    glats = 2 * dc.msl.attrs['GRIB_N']
    glons = 2 * glats
    
    dlon = 360./glons
    
    npoints = (maxlon - minlon)/dlon + 1
    
    x = np.linspace(minlon,maxlon,npoints)
    
    values, index = np.unique(dc.latitude,return_index=True)
    
    new_data = []
    for val in values:
        row = dc.where(dc.latitude==val, drop=True)
        da = row.assign_coords(values=row.longitude).interp(values=x,kwargs={'fill_value': 'extrapolate'}, method='cubic')
        da = da.drop('longitude').rename({'values':'longitude'}).drop('latitude')
        new_data.append(da)
    #concat
    res = xr.concat(new_data, dim='latitude').transpose('time', 'latitude', 'longitude')
        
    res = res.assign_coords(latitude = values)
            
    logger.info('regriding done')
    
    
    return res
    
    
def regrid(ds):
    
    logger.info('regriding meteo')

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
    
    logger.info('regriding done')
    
    
    return data



class meteo:
   
    def __init__(self, meteo_files=None, engine=None, url=None, **kwargs):
        
        """Read meteo data from variable sources.
 
        Parameters
        ----------
        mfiles : str
            list of files
        engine : str
            Name of xarray backend to be used
        url : str
            url for an online server (erdapp, etc.)
        
        Returns
        -------
        retrieved : xarray DataSet

        """
                     
        if meteo_files:                    
            if engine == 'cfgrib' :
                self.Dataset = cfgrib(meteo_files, **kwargs)
            elif engine == 'pynio' :
                self.Dataset = pynio(meteo_files, **kwargs)
            elif engine == 'netcdf' :
                self.Dataset = netcdf(meteo_files, **kwargs)
            else:
                if 'grib' in meteo_files[0].split('.')[-1]:
                    self.Dataset = cfgrib(meteo_files, **kwargs)
                elif 'nc' in meteo_files[0].split('.')[-1]:  
                    self.Dataset = netcdf(meteo_files, **kwargs)
                           
        else:        
            self.Dataset = from_url(**kwargs)
                        
        

def cfgrib(filenames=None, minlon=None, maxlon=None, minlat=None, maxlat=None, start_date=None, end_date=None, time_frame=None, irange=[0,-1,1], combine=False, **kwargs):

    backend_kwargs = kwargs.get('backend_kwargs', {'indexpath':''})
    xr_kwargs = kwargs.get('xr_kwargs', {'concat_dim':'step'})

    try:
        start_date = pd.to_datetime(start_date)
    except:
        pass

    if time_frame:
        try:
            end_date = start_date + pd.to_timedelta(time_frame)
        except:
            pass            
    else:
        try:
            end_date = pd.to_datetime(end_date)
        except:
            pass

    ft1, ft2, dft = irange
        
    ts = pd.to_datetime(start_date)
    te = pd.to_datetime(end_date)     
    
    #--------------------------------------------------------------------- 
    logger.info('extracting meteo')
    #---------------------------------------------------------------------      

    data = xr.open_mfdataset(filenames, combine='nested', engine='cfgrib', backend_kwargs=backend_kwargs, **xr_kwargs)    

    data = data.squeeze(drop=True)
    #        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing> not efficient for output  
    
    time_coord = data.msl.dims[0]

    if time_coord != 'time':
        if 'time' in data.coords.keys():
            data = data.rename({'time':'rtime'})
            data = data.rename({time_coord:'time'})
            data = data.assign_coords(time=data.valid_time)
        
    if combine :
        mask = data.time.to_pandas().duplicated('last').values 
        msl = data.msl[~mask]
        u10 = data.u10[~mask]
        v10 = data.v10[~mask]
        data = xr.merge([msl,u10,v10])
           
        
    if not minlon : minlon = data.longitude.data.min()
    if not maxlon : maxlon = data.longitude.data.max()
    if not minlat : minlat = data.latitude.data.min()
    if not maxlat : maxlat = data.latitude.data.max()


    if minlon < data.longitude.data.min() : minlon = minlon + 360.

    if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.

    if minlon > data.longitude.data.max() : minlon = minlon - 360.

    if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.
                    
    
    if not ts : ts = data.time.data[ft1]
    if not te : te = data.time.data[ft2]
        

    if ts < data.time.data.min() :
        logger.warning('coverage between {} and {} \n'.format(pd.to_datetime(data.valid_time.min().values).strftime('%Y.%m.%d %H:%M:%S'),pd.to_datetime(data.valid_time.max().values).strftime('%Y.%m.%d %H:%M:%S')))
        logger.error('time frame does not match source range\n')
        sys.exit(1)
  
    if te > data.time.data.max() :
        logger.warning('coverage between {} and {} \n'.format(pd.to_datetime(data.valid_time.min().values).strftime('%Y.%m.%d %H:%M:%S'),pd.to_datetime(data.valid_time.max().values).strftime('%Y.%m.%d %H:%M:%S')))
        logger.error('time frame does not match source range\n')
        sys.exit(1)


    if len(data.longitude.shape) == 2:
        d1 = data.where(data.longitude>minlon,drop=True)
        d2 = d1.where(d1.longitude<maxlon,drop=True)
        d3 = d2.where(d2.latitude>minlat,drop=True)
        d4 = d3.where(d3.latitude<maxlat,drop=True)        
    
        data = regrid(d4)
    
    if data.msl.attrs['GRIB_gridType'] == 'reduced_gg':
        d1 = data.where(data.longitude>minlon,drop=True)
        d2 = d1.where(d1.longitude<maxlon,drop=True)
        d3 = d2.where(d2.latitude>minlat,drop=True)
        d4 = d3.where(d3.latitude<maxlat,drop=True)        
    
        data = reduced_gg(d4)
    

    tslice=slice(ts, te, dft)

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
      
    # TODO
    #        if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
    #            c = np.sign(np.mean(minlon, maxlon))    
    #            tot.longitude = tot.longitude + c*360.


    #--------------------------------------------------------------------- 
    logger.info('meteo done\n')
    #--------------------------------------------------------------------- 
    
    return tot


    
def pynio(filenames=None, minlon=None, maxlon=None, minlat=None, maxlat=None, start_date=None, end_date=None, time_frame=None, irange=[0,-1,1], combine=False, **kwargs):
    
    backend_kwargs = kwargs.get('backend_kwargs', {})
    xr_kwargs = kwargs.get('xr_kwargs', {'concat_dim':'step'})
    if 'preprocess' in xr_kwargs.keys():
        xr_kwargs['preprocess'] = fix_my_data        

    try:
        start_date = pd.to_datetime(start_date)
    except:
        pass

    if time_frame:
        try:
            end_date = start_date + pd.to_timedelta(time_frame)
        except:
            pass            
    else:
        try:
            end_date = pd.to_datetime(end_date)
        except:
            pass

    ft1, ft2, dft = irange


    ts = pd.to_datetime(start_date)
    te = pd.to_datetime(end_date)     

    #--------------------------------------------------------------------- 
    logger.info('extracting meteo')
    #---------------------------------------------------------------------      

    data = xr.open_mfdataset(filenames, combine='nested', engine='pynio', backend_kwargs=backend_kwargs, **xr_kwargs)    

    data = data.squeeze(drop=True)
        
    time_coord = [x for x in data.coords if 'time' in data[x].long_name.lower() ]
    lon_coord = [x for x in data.coords if 'longitude' in data[x].long_name.lower() ]
    lat_coord = [x for x in data.coords if 'latitude' in data[x].long_name.lower() ]
    msl_ = [x for x in data.variables if 'pressure' in data[x].long_name.lower() ]
    u10_ = [x for x in data.variables if ('u wind' in data[x].long_name.lower()) | ('u-component' in data[x].long_name.lower())]
    v10_ = [x for x in data.variables if ('v wind' in data[x].long_name.lower()) | ('v-component' in data[x].long_name.lower())]
 
    data = data.rename({msl_[0]:'msl',u10_[0]:'u10',v10_[0]:'v10',lon_coord[0]:'longitude',lat_coord[0]:'latitude'})

    if time_coord:
        name = time_coord[0]
        if 'forecast' in name:               
            tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(data.coords[name].values)
            data = data.assign_coords(name=tts)
        data = data.rename({name:'time'})
    else:
        tts = pd.to_datetime(data.msl.attrs['initial_time'], format='%m/%d/%Y (%H:%M)') + pd.to_timedelta(data.coords['step'].values, unit='H')
        data = data.rename({'step':'time'})
        data = data.assign_coords(time=tts)
    
#    if combine : TODO
#        mask = data.time.to_pandas().duplicated('last').values 
#        msl = data.msl[~mask]
#        u10 = data.u10[~mask]
#        v10 = data.v10[~mask]
#        data = xr.merge([msl,u10,v10])
    

    #        data = data.sortby('latitude', ascending=True)   # make sure that latitude is increasing      
                
    if not minlon : minlon = data.longitude.data.min()
    if not maxlon : maxlon = data.longitude.data.max()
    if not minlat : minlat = data.latitude.data.min()
    if not maxlat : maxlat = data.latitude.data.max()


    if minlon < data.longitude.data.min() : minlon = minlon + 360.

    if maxlon < data.longitude.data.min() : maxlon = maxlon + 360.

    if minlon > data.longitude.data.max() : minlon = minlon - 360.

    if maxlon > data.longitude.data.max() : maxlon = maxlon - 360.

   
    
    if not ts : ts = data.time.data[ft1]
    if not te : te = data.time.data[ft2]


    if ts < data.time.data.min() :
        logger.warning('coverage between {} and {} \n'.format(data.time.data.min(),data.time.data.max()))
        logger.error('time frame not available\n')
        sys.exit(1)
  
    if te > data.time.data.max() :
        logger.error('time frame not available\n')
        logger.warning('coverage between {} and {} \n'.format(data.time.data.min(),data.time.data.max()))
        sys.exit(1)

    if len(data.longitude.shape) == 2:
        d1 = data.where(data.longitude>minlon,drop=True)
        d2 = d1.where(d1.longitude<maxlon,drop=True)
        d3 = d2.where(d2.latitude>minlat,drop=True)
        d4 = d3.where(d3.latitude<maxlat,drop=True)        
    
        data = regrid(d4)

    tslice=slice(ts, te, dft)

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

# TODO    
#    if np.abs(np.mean([minlon,maxlon]) - np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)])) > 300. :
#        c = np.sign(np.mean([kwargs.get('minlon', None), kwargs.get('maxlon', None)]))    
#        tot.longitude = tot.longitude + c*360.

    #--------------------------------------------------------------------- 
    logger.info('meteo done\n')
    #--------------------------------------------------------------------- 
    
    return tot



def from_url(url = None, minlon=None, maxlon=None, minlat=None, maxlat=None, start_date=None, end_date=None, time_frame=None, **kwargs):
        
    try:
        start_date = pd.to_datetime(start_date)
    except:
        pass

    if time_frame:
        try:
            end_date = start_date + pd.to_timedelta(time_frame)
        except:
            pass            
    else:
        try:
            end_date = pd.to_datetime(end_date)
        except:
            pass
    
    ts = pd.to_datetime(start_date)
    te = pd.to_datetime(end_date)     

    
    if not url:
        for link in erddap():
            try:
                data = xr.open_dataset(link)
                url = link
                break
            except:
                continue
                 
    #--------------------------------------------------------------------- 
    logger.info('extracting meteo from {}.html\n'.format(url))
    #---------------------------------------------------------------------      
        
    data = xr.open_dataset(url)  
    
    try:
        data = data.rename({'lon':'longitude','lat':'latitude'})
    except:
        pass  
 
    lon0 = minlon + 360. if minlon < data.longitude.min() else minlon
    lon1 = maxlon + 360. if maxlon < data.longitude.min() else maxlon

    lon0 = lon0 - 360. if lon0 > data.longitude.max() else lon0
    lon1 = lon1 - 360. if lon1 > data.longitude.max() else lon1
    
    if ts < data.time.min().values :
      ld = pd.to_datetime(data.time.min().values).strftime(format='%Y-%m-%d %H:%M:%S')
      hd = pd.to_datetime(data.time.max().values).strftime(format='%Y-%m-%d %H:%M:%S')
      logger.error('time frame not available\n')
      logger.warning('coverage between {} and {} \n'.format(ld,hd))
      sys.exit(1)
  
    if te > data.time.max().values :
      logger.error('time frame not available\n')
      logger.warning('coverage between {} and {} \n'.format(ld,hd))
      sys.exit(1)

    tslice=slice(ts, te)    

    i0=np.abs(data.longitude.data-lon0).argmin()
    i1=np.abs(data.longitude.data-lon1).argmin()


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
      
    if np.abs(np.mean(tot.longitude) - np.mean([minlon, maxlon])) > 300. :
      c = np.sign(np.mean([minlon, maxlon]))    
      tot['longitude'] = tot['longitude'] + c*360.



    tot = tot.rename({'prmslmsl':'msl','ugrd10m':'u10','vgrd10m':'v10'})
     
    #--------------------------------------------------------------------- 
    logger.info('meteo done\n')
    #--------------------------------------------------------------------- 

    return tot
        


def netcdf(filename=None, **kwargs):

        
    #--------------------------------------------------------------------- 
    logger.info('extracting meteo\n')
    #---------------------------------------------------------------------      
        
    return xr.open_mfdataset(filenames, combine='nested')
         
    #--------------------------------------------------------------------- 
    logger.info('meteo done\n')
    #--------------------------------------------------------------------- 

    

def to_output(dataset=None,solver=None, **kwargs):
                        
    model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
    s = getattr(model,solver) # get solver class
                
    s.to_force(dataset,vars=['msl','u10','v10'], **kwargs)
        
       

    
    