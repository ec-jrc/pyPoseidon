"""
Dem module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import pyresample
import xarray as xr
import sys
import importlib
from pyPoseidon.utils.fix import fix
import logging


logger = logging.getLogger('pyPoseidon')

class dem:
    def __init__(self, dem_source=None, **kwargs):
        if not dem_source : dem_source = 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus'
        self.Dataset = dem_(source = dem_source, **kwargs)
        
        coastline = kwargs.get('coastlines', None) 
        if coastline is not None :
            if kwargs.get('adjust_dem', True) : self.adjust(coastline)   
     
    def adjust(self,coastline,**kwargs):
         
        self.Dataset = fix(self.Dataset,coastline,**kwargs)
  
    
      
def dem_(source=None, lon_min=-180, lon_max=180, lat_min=-90, lat_max=90, **kwargs):
    
    ncores = kwargs.get('ncores', 1)       
 
    #--------------------------------------------------------------------- 
    logger.info('extracting dem from {}\n'.format(source))
    #---------------------------------------------------------------------      
  
    data = xr.open_dataset(source) 
    
    #rename vars,coords
    var = [keys for keys in data.data_vars]
    coords = [keys for keys in data.coords]
    lat = [x for x in coords if 'lat' in  x]  
    lon = [x for x in coords if 'lon' in  x]
    data = data.rename({var[0] : 'elevation', lat[0] : 'latitude', lon[0] : 'longitude'})
    
    #recenter the window 
    
    lon0 = lon_min + 360. if lon_min < -180 else lon_min
    lon1 = lon_max + 360. if lon_max < -180 else lon_max

    lon0 = lon0 - 360. if lon0 > 180 else lon0
    lon1 = lon1 - 360. if lon1 > 180 else lon1
    
#   TODO check this for regional files
    if (lon_min < data.longitude.min()) or (lon_max > data.longitude.max()): 
        logger.warning('Lon must be within {} and {}'.format(data.longitude.min().values,data.longitude.max().values))
        logger.warning('compensating if global dataset available')
        
#        sys.exit()
    if (lat_min < data.latitude.min()) or (lat_max > data.latitude.max()): 
        logger.warning('Lat must be within {} and {}'.format(data.latitude.min().values,data.latitude.max().values))
        logger.warning('compensating if global dataset available')
#        sys.exit()

    #get idx
    i0=np.abs(data.longitude.data-lon0).argmin()
    i1=np.abs(data.longitude.data-lon1).argmin()


    j0=np.abs(data.latitude.data-lat_min).argmin()
    j1=np.abs(data.latitude.data-lat_max).argmin()


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

        p1 = (
          data.elevation
          .isel(longitude=slice(lon_0,data.longitude.size),latitude=slice(lat_0,lat_1))
          )

        p1 = p1.assign_coords({'longitude':p1.longitude.values - 360.})


        p2 = (
          data.elevation
          .isel(longitude=slice(0,lon_1),latitude=slice(lat_0,lat_1))
          )

        dem = xr.concat([p1,p2],dim='longitude')
  
    else:            

        dem = (
            data.elevation
            .isel(longitude=slice(lon_0,lon_1),latitude=slice(lat_0,lat_1))
            )
    
    
    if np.abs(np.mean(dem.longitude) - np.mean([lon_min, lon_max])) > 170. :
        c = np.sign(np.mean([lon_min, lon_max]))    
        dem['longitude'] = dem['longitude'] + c*360.
        


    if 'grid_x' in kwargs.keys():

        #--------------------------------------------------------------------- 
        logger.info('.. interpolating on grid ..\n')
        #---------------------------------------------------------------------      

        grid_x = kwargs.get('grid_x', None)
        grid_y = kwargs.get('grid_y', None)
        # resample on the given grid
        xx,yy = np.meshgrid(dem.longitude ,dem.latitude)   #original grid         

        # Translate for pyresample
        if xx.mean() < 0 and xx.min() < -180. :
            xx = xx + 180.
            gx = grid_x + 180.
        elif xx.mean() > 0 and xx.max() > 180. :
            xx = xx - 180.
            gx = grid_x - 180.
        else:
            gx = grid_x
        

        orig = pyresample.geometry.SwathDefinition(lons=xx,lats=yy) # original points
        targ = pyresample.geometry.SwathDefinition(lons=gx,lats=grid_y) # target grid

        
        wet = kwargs.get('wet_only', False)
        if wet:
            #mask positive bathymetry 
            vals = np.ma.masked_array(dem,dem.values > 0 )
        else:
            vals = dem.values


        # with nearest using only the water values        

        #    itopo = pyresample.kd_tree.resample_nearest(orig,dem.values,targ,radius_of_influence=50000,fill_value=np.nan,nprocs=ncores)

        grid_con = pyresample.image.ImageContainerNearest(vals, orig, radius_of_influence=50000,fill_value=np.nan)#,nprocs=ncores)

        area_con = grid_con.resample(targ)

        itopo = area_con.image_data

        if len(grid_x.shape) > 1:         
            idem = xr.Dataset({'ival': (['k', 'l'],  itopo), 
                           'ilons': (['k', 'l'], grid_x),   
                           'ilats': (['k', 'l'], grid_y)})#, 
#                           coords={'ilon': ('ilon', grid_x[0,:]),   
#                                   'ilat': ('ilat', grid_y[:,0])}) 
                                           
        elif len(grid_x.shape) == 1:
            idem = xr.Dataset({'ival': (['k'],  itopo), 
                    'ilons': (['k'], grid_x),   
                    'ilats': (['k'], grid_y)}
                         )

        #--------------------------------------------------------------------- 
        logger.info('dem done\n')
        #--------------------------------------------------------------------- 
         

        return xr.merge([dem,idem])
        
    else:
        
        return xr.merge([dem])
 

 
def to_output(dataset=None,solver=None, **kwargs):
                
    model=importlib.import_module('pyPoseidon.model') #load pyPoseidon model class

    s = getattr(model,solver) # get solver class
    
    s.to_dep(dataset, **kwargs)
       
    