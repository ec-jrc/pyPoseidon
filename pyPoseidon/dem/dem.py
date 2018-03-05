import numpy as np
import netCDF4
import scipy.interpolate
import pyresample
from pyPoseidon.utils.bfs import *
from mpl_toolkits.basemap import Basemap
from shapely import geometry, ops
import matplotlib.path as mpltPath
import geopandas as gp
import xarray as xr
import sys

class dem:
    impl=None
    def __init__(self, **kwargs):
        dem = kwargs.get('dem', None)
        if dem == 'gebco' :
                self.impl = gebco(**kwargs)    
        elif dem == 'emodnet' :
                self.impl = emodnet(**kwargs)    
        else:
            self.impl = erdap(**kwargs)


class emodnet(dem):

    def __init__(self,**kwargs):
    
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None)       
      properties = kwargs.get('properties', {})    
        
      url = kwargs.get('dpath', None)      
    # open NetCDF data in 
      data = xr.open_dataset(url)    
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting dem from {}\n'.format(url))
      sys.stdout.flush()
      #---------------------------------------------------------------------      
      
      lon=data.longitude.data
      lat=data.latitude.data     
      
      if (minlon < lon.min()) or (maxlon > lon.max()): print 'Lon must be within {} and {}'.format(lon.min(),lon.max())
      if (minlat < lat.min()) or (maxlat > lat.max()): print 'Lat must be within {} and {}'.format(lat.min(),lat.max())

      i1=np.abs(lon-minlon).argmin()
      if i1 > 0: i1=i1-1
      i2=np.abs(lon-maxlon).argmin()
      if i2 < lon.shape : i2=i2+1

      j1=np.abs(lat-minlat).argmin()
      if j1 > 0: j1=j1-1
      j2=np.abs(lat-maxlat).argmin()
      if j2 < lat.shape: j2=j2+1

      lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
      topo = data.depth[j1:j2,i1:i2]
        
      self.val = topo
      self.dlons = lons 
      self.dlats = lats
         
      if 'grid_x' in kwargs.keys():
       grid_x = kwargs.get('grid_x', None)
       grid_y = kwargs.get('grid_y', None)
    # resample on the given grid
       ilon=lons[0,:]
       ilat=lats[:,0]
              
       orig = pyresample.geometry.SwathDefinition(lons=self.dlons,lats=self.dlats) # original points
       targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
       itopo = pyresample.kd_tree.resample_nearest(orig,self.val.data,targ,radius_of_influence=50000,fill_value=999999)

       self.ival = itopo
       self.ilons = grid_x
       self.ilats = grid_y
       
       
       
       #--------------------------------------------------------------------- 
       sys.stdout.flush()
       sys.stdout.write('\n')
       sys.stdout.write('dem done\n')
       sys.stdout.flush()
       #--------------------------------------------------------------------- 
       

class erdap(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
           
      if minlon < -180 : minlon = minlon + 360.    
      if maxlon < -180 : maxlon = maxlon + 360.
      if minlon >  180 : minlon = minlon - 360.    
      if maxlon >  180 : maxlon = maxlon - 360.
      
      
      
      url = kwargs.get('url', 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus')
            
    #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting dem from {}\n'.format(url))
      sys.stdout.flush()
    #---------------------------------------------------------------------      
          
            
      data = xr.open_dataset(url)    
      
      i0=np.abs(data.longitude.data-minlon).argmin()
      i1=np.abs(data.longitude.data-maxlon).argmin()

      
      j0=np.abs(data.latitude.data-minlat).argmin()
      j1=np.abs(data.latitude.data-maxlat).argmin()
      
      
      if i0 > i1 :

          sh = (
              data.z
              .isel(longitude=slice(i0,data.longitude.size),latitude=slice(j0,j1))
              )
          sh.longitude.values = sh.longitude.values -360.

          sh1 = (
              data.z
              .isel(longitude=slice(0,i1),latitude=slice(j0,j1))
              )
              
          dem = xr.concat([sh,sh1],dim='longitude')
          
      else:            
      
          dem = (
              data.z
              .isel(longitude=slice(i0,i1),latitude=slice(j0,j1))
                 )
      
      self.val = dem
      
      
      
      if 'grid_x' in kwargs.keys():
         grid_x = kwargs.get('grid_x', None)
         grid_y = kwargs.get('grid_y', None)
      # resample on the given grid
      
         xx, yy = np.meshgrid(self.val.longitude.data,self.val.latitude.data)
              
         orig = pyresample.geometry.SwathDefinition(lons=xx,lats=yy) # original points
         targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
         itopo = pyresample.kd_tree.resample_nearest(orig,self.val.data,targ,radius_of_influence=50000,fill_value=999999)

         self.ival = itopo
         self.ilons = grid_x
         self.ilats = grid_y
      
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('dem done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
class gebco(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
           
      if minlon < -180: minlon = minlon + 360.
      
      if maxlon < -180: maxlon = maxlon + 360.
    
      url = kwargs.get('dpath', None)
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting dem from {}\n'.format(url))
      sys.stdout.flush()
      #---------------------------------------------------------------------      
      
            
      data = xr.open_dataset(url)    
      
      i0=np.abs(data.lon.data-minlon).argmin()
      i1=np.abs(data.lon.data-maxlon).argmin()

      
      j0=np.abs(data.lat.data-minlat).argmin()
      j1=np.abs(data.lat.data-maxlat).argmin()
      
      dem = (
          data[data.data_vars.keys()[0]]
          .isel(lon=slice(i0,i1),lat=slice(j0,j1))
          )
      
      self.val = dem
      
      
      if 'grid_x' in kwargs.keys():
         grid_x = kwargs.get('grid_x', None)
         grid_y = kwargs.get('grid_y', None)
      # resample on the given grid
      
         xx, yy = np.meshgrid(self.val.lon.data,self.val.lat.data)
              
         orig = pyresample.geometry.SwathDefinition(lons=xx,lats=yy) # original points
         targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
         itopo = pyresample.kd_tree.resample_nearest(orig,self.val.data,targ,radius_of_influence=50000,fill_value=999999)

         self.ival = itopo
         self.ilons = grid_x
         self.ilats = grid_y
      
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('dem done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
     

    