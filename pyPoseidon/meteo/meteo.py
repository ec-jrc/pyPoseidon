# -*- coding: utf-8 -*- 
"""
Meteo module. Pre-processing the weather forcing component.

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
   
    def __init__(self, mfiles=None, engine=None, url=None, **kwargs):
        
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
             
        if not url : url = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/NCEP_Global_Best'
        
        if mfiles:                    
            if engine == 'cfgrib' :
                self.uvp = cfgrib(mfiles, **kwargs)
            elif engine == 'pynio' :
                self.uvp = pynio(mfiles, **kwargs)
            elif engine == 'netcdf' :
                self.uvp = netcdf(mfiles, **kwargs)         
        else:        
            self.uvp = from_url(url=url, **kwargs)
                        
        

def cfgrib(filenames=None, minlon=None, maxlon=None, minlat=None, maxlat=None, range=None, ft1=0, ft2=-1, dft=1, combine=False, **kwargs):

    backend_kwargs = kwargs.get('backend_kwargs', {'indexpath':''})
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

    data = xr.open_mfdataset(filenames, engine='cfgrib', backend_kwargs=backend_kwargs, **xr_kwargs)    

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
                    
    
    if not ts : ts = data.time.data.min()
    if not te : te = data.time.data.max()

    ts = data.time.data[ft1]

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


    
def pynio(filenames=None, minlon=None, maxlon=None, minlat=None, maxlat=None, range=None, ft1=0, ft2=-1, dft=1, combine=False, **kwargs):
    
    backend_kwargs = kwargs.get('backend_kwargs', {})
    xr_kwargs = kwargs.get('xr_kwargs', {'concat_dim':'step'})
    if 'preprocess' in xr_kwargs.keys():
        xr_kwargs['preprocess'] = fix_my_data        


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

    data = xr.open_mfdataset(filenames, engine='pynio', backend_kwargs=backend_kwargs, **xr_kwargs)    

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

   
    
    if not ts : ts = data.time.data.min()
    if not te : te = data.time.data.max()

    ts = data.time.data[ft1]

    te = data.time.data[ft2]


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



def from_url(url = None, minlon=None, maxlon=None, minlat=None, maxlat=None, **kwargs):
                            
    ts = kwargs.get('start_date', None)
    te = kwargs.get('end_date', None)

    ts = pd.to_datetime(ts)
    te = pd.to_datetime(te)     

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
      logger.error('time frame not available\n')
      logger.warning('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
      sys.exit(1)
  
    if te > data.time.max().values :
      logger.error('time frame not available\n')
      logger.warning('coverage between {} and {} \n'.format(data.time.min(),data.time.max()))
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
        
    return xr.open_mfdataset(filenames)
         
    #--------------------------------------------------------------------- 
    logger.info('meteo done\n')
    #--------------------------------------------------------------------- 

    

def to_output(dataset=None,solver=None, **kwargs):
                        
    model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class
                
    s = getattr(model,solver) # get solver class
                
    s.to_force(dataset,vars=['msl','u10','v10'], **kwargs)
        
    
def read_d3d_meteo(filename=None, name=None):
    
    df = pd.read_csv(filename,header=0, names=['data'], index_col=None, low_memory=False)
    
    tlines = df[df.data.str.contains('TIME')].index # rows which start with TIME
    
    # get attrs
    d1 = df.loc[0:tlines[0]-1,'data'].str.split('=', 2, expand=True)
    d1.columns=['key','value'] # assign column names
    d1.key = d1.key.str.strip() # cleanup spaces
    d1.value = d1.value.str.strip()
    attrs = dict(zip(d1.key, d1.value)) # create dict
    for key in ['n_cols','n_rows','n_quantity']: # str -> int
        attrs[key] = int(attrs[key])
        
    for key in ['x_llcenter','dx','y_llcenter','dy','NODATA_value']:
        attrs[key] = float(attrs[key])
    
    # get time reference
    d2 = df.loc[tlines, 'data'].str.split('=', 2, expand=True)
    d2 = d2.drop(d2.columns[0], axis=1)
    d2.columns=['data']
    d2 = d2.loc[:,'data'].str.split(' ', 4, expand=True)
    d2 = d2.drop(d2.columns[[0,2,3]], axis=1)
    d2.columns = ['hours','time0']
    d2.hours = d2.hours.apply(pd.to_numeric)
    d2.time0 = pd.to_datetime(d2.time0.values)
    d2 = d2.reset_index(drop=True)
    #create timestamps
    time=[]
    for i in range(d2.shape[0]):
        time.append(d2.time0[0] + pd.DateOffset(hours=int(d2.loc[i,'hours'])))
    d2['time']=time
    
    #get the float numbers
    d3 = df.drop(np.arange(0,tlines[0]))
    d3 = d3.drop(tlines)
    
#    data = []
#    for i in range(d3.values.shape[0]):
#        row = d3.values[i][0].split(' ')
#        row = [np.float(x) for x in row]
#        data.append(row)
#    data = np.array(data) # make array
    
    data = d3[d3.columns[0]].str.split(' ', attrs['n_cols'], expand=True).to_numpy().astype(float)
    
    data = data.reshape(d2.shape[0],attrs['n_rows'], attrs['n_cols']) # reshape
    
    #define lat/lon
    lon = [attrs['x_llcenter'] + attrs['dx'] * i for i in np.arange(attrs['n_cols'])]
    lat = [attrs['y_llcenter'] + attrs['dy'] * i for i in np.arange(attrs['n_rows'])]
    
    #create an xarray
    da = xr.DataArray(data, dims=['time','latitude','longitude'],
                         coords={'time': d2.time, 'latitude':lat, 'longitude':lon}, name=name)
    
    da.attrs = attrs
                    
    return da
    
    

    
    