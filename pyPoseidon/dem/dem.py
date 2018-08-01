import numpy as np
import pyresample
import xarray as xr
import sys
from pyPoseidon.utils.fix import *

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
      ncores = kwargs.get('ncores', 1)       
      properties = kwargs.get('properties', {})    
        
      url = kwargs.get('dpath', None)      
    # open NetCDF data in 
      data = xr.open_mfdataset(url)    
      
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
              
      self.Dataset = xr.Dataset({'val': (['dlat', 'dlon'],  topo), 
                            'dlons': (['i', 'j'], lons),   
                            'dlats': (['i', 'j'], lats)}, 
                            coords={'dlon': ('dlon', lons[0,:]),   
                                    'dlat': ('dlat', lats[:,0])})         
      
         
      if 'grid_x' in kwargs.keys():
          #--------------------------------------------------------------------- 
          sys.stdout.flush()
          sys.stdout.write('\n')
          sys.stdout.write('.. interpolating on grid ..\n')
          sys.stdout.flush()
          #---------------------------------------------------------------------      
          
          
          grid_x = kwargs.get('grid_x', None)
          grid_y = kwargs.get('grid_y', None)
           # resample on the given grid
              
          orig = pyresample.geometry.SwathDefinition(lons=lons,lats=lats) # original points
          targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
 #         itopo = pyresample.kd_tree.resample_nearest(orig,topo.values,targ,radius_of_influence=50000,fill_value=np.nan, nprocs=ncores)

          grid_con = pyresample.image.ImageContainerNearest(topo.values, orig, radius_of_influence=50000,fill_value=np.nan,nprocs=ncores)
 
          area_con = grid_con.resample(targ)

          itopo = area_con.image_data

       
          if len(grid_x.shape) > 1: 
                 
              dem = xr.Dataset({'ival': (['ilat', 'ilon'],  itopo), 
                               'ilons': (['k', 'l'], grid_x),   
                               'ilats': (['k', 'l'], grid_y)}, 
                               coords={'ilon': ('ilon', grid_x[0,:]),   
                                       'ilat': ('ilat', grid_y[:,0])})         
          elif len(grid_x.shape) == 1:
             
             dem = xr.Dataset({'ival': (['k'],  itopo), 
                        'ilons': (['k'], grid_x),   
                        'ilats': (['k'], grid_y)}
                             )
      
          self.Dataset = xr.merge([self.Dataset,dem])
       
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('dem done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
       
    def adjust(self,shpfile,**kwargs):
         
         wmask, cg = fix(self,shpfile,**kwargs)
         bmatch(self,wmask,**kwargs)



class erdap(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
      ncores = kwargs.get('ncores', 1)       
      
           
      if minlon < -180 : minlon = minlon + 360.    
      if maxlon < -180 : maxlon = maxlon + 360.
      if minlon >  180 : minlon = minlon - 360.    
      if maxlon >  180 : maxlon = maxlon - 360.
      
      
      
      url = kwargs.get('dem_url', 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus')
            
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
            
      self.Dataset = dem       
      
      
      if 'grid_x' in kwargs.keys():
         #--------------------------------------------------------------------- 
         sys.stdout.flush()
         sys.stdout.write('\n')
         sys.stdout.write('.. interpolating on grid ..\n')
         sys.stdout.flush()
         #---------------------------------------------------------------------      
          
          
         grid_x = kwargs.get('grid_x', None)
         grid_y = kwargs.get('grid_y', None)
      # resample on the given grid
      
         xx, yy = np.meshgrid(dem.longitude.data,dem.latitude.data)
              
         orig = pyresample.geometry.SwathDefinition(lons=xx,lats=yy) # original points
         targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
#         itopo = pyresample.kd_tree.resample_nearest(orig,dem.values,targ,radius_of_influence=50000,fill_value=np.nan, nprocs=ncores)
 
         grid_con = pyresample.image.ImageContainerNearest(dem.values, orig, radius_of_influence=50000,fill_value=np.nan,nprocs=ncores)

         area_con = grid_con.resample(targ)

         itopo = area_con.image_data


         if len(grid_x.shape) > 1: 
             idem = xr.Dataset({'ival': (['ilat', 'ilon'],  itopo), 
                               'ilons': (['k', 'l'], grid_x),   
                               'ilats': (['k', 'l'], grid_y)}, 
                               coords={'ilon': ('ilon', grid_x[0,:]),   
                                       'ilat': ('ilat', grid_y[:,0])})         
      
         elif len(grid_x.shape) == 1:
            idem = xr.Dataset({'ival': (['k'],  itopo), 
                              'ilons': (['k'], grid_x),   
                              'ilats': (['k'], grid_y)} 
                              )         
            
      
         self.Dataset = xr.merge([self.Dataset,idem])
               
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('dem done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
      
    def adjust(self,shpfile,**kwargs):
         
         wmask, cg = fix(self,shpfile,**kwargs)
         bmatch(self,wmask,**kwargs)
         
      
class gebco(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
      ncores = kwargs.get('ncores', 1)       
      
                 
      if minlon < -180: minlon = minlon + 360.
      
      if maxlon < -180: maxlon = maxlon + 360.
    
      url = kwargs.get('dpath', None)
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('extracting dem from {}\n'.format(url))
      sys.stdout.flush()
      #---------------------------------------------------------------------      
      
            
      data = xr.open_mfdataset(url)    
      
      i0=np.abs(data.lon.data-minlon).argmin()
      i1=np.abs(data.lon.data-maxlon).argmin()

      
      j0=np.abs(data.lat.data-minlat).argmin()
      j1=np.abs(data.lat.data-maxlat).argmin()
      
      dem = (
          data[data.data_vars.keys()[0]]
          .isel(lon=slice(i0,i1),lat=slice(j0,j1))
          )
      
      xx,yy = np.meshgrid(dem.lon,dem.lat)
      
      
      self.Dataset = xr.Dataset({'val': (['dlat', 'dlon'],  dem), 
                            'dlons': (['i', 'j'], xx),   
                            'dlats': (['i', 'j'], yy)}, 
                            coords={'dlon': ('dlon', xx[0,:]),   
                                    'dlat': ('dlat', yy[:,0])})         
      
      
      if 'grid_x' in kwargs.keys():
          
         #--------------------------------------------------------------------- 
         sys.stdout.flush()
         sys.stdout.write('\n')
         sys.stdout.write('.. interpolating on grid ..\n')
         sys.stdout.flush()
         #---------------------------------------------------------------------      
          
         grid_x = kwargs.get('grid_x', None)
         grid_y = kwargs.get('grid_y', None)
      # resample on the given grid
                    
         orig = pyresample.geometry.SwathDefinition(lons=xx,lats=yy) # original points
         targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
         # with nearest using only the water values        
       
     #    itopo = pyresample.kd_tree.resample_nearest(orig,dem.values,targ,radius_of_influence=50000,fill_value=np.nan,nprocs=ncores)
         
         grid_con = pyresample.image.ImageContainerNearest(dem.values, orig, radius_of_influence=50000,fill_value=np.nan,nprocs=ncores)

         area_con = grid_con.resample(targ)

         itopo = area_con.image_data
        
         if len(grid_x.shape) > 1:         
             dem = xr.Dataset({'ival': (['ilat', 'ilon'],  itopo), 
                               'ilons': (['k', 'l'], grid_x),   
                               'ilats': (['k', 'l'], grid_y)}, 
                               coords={'ilon': ('ilon', grid_x[0,:]),   
                                       'ilat': ('ilat', grid_y[:,0])})         
         elif len(grid_x.shape) == 1:
             dem = xr.Dataset({'ival': (['k'],  itopo), 
                        'ilons': (['k'], grid_x),   
                        'ilats': (['k'], grid_y)}
                             )
                         
       
         self.Dataset = xr.merge([self.Dataset,dem])
         
      
      #--------------------------------------------------------------------- 
      sys.stdout.flush()
      sys.stdout.write('\n')
      sys.stdout.write('dem done\n')
      sys.stdout.flush()
      #--------------------------------------------------------------------- 
      
     
    def adjust(self,shpfile,**kwargs):
         
         wmask, cg = fix(self,shpfile,**kwargs)
         bmatch(self,wmask,**kwargs)
         
    