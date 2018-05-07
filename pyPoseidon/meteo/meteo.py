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

            grib_release(gid)

            return name,dat,lon,lat


class meteo:
    impl=None
    def __init__(self,**kwargs):
        msource = kwargs.get('meteo', None)
        if msource == 'ecmwf_oper' :
            self.impl = ecmwf_oper(**kwargs)
        elif msource == 'hnms_oper' :
            self.impl = hnms_oper(**kwargs)
        elif msource == 'am_oper' :
            self.impl = am_oper(**kwargs)
        elif msource == 'generic' :
            self.impl = generic(**kwargs)
        elif msource == 'gfs_oper' :
            self.impl = gfs_oper(**kwargs)
        
        else:
            self.impl = gfs_erdap(**kwargs)
                
    
class ecmwf_oper(meteo):   
        
    def __init__(self,**kwargs):
    
      filenames = kwargs.get('mpaths', {})
      ft1 = kwargs.get('ft1', 0)
      ft2 = kwargs.get('ft2', 11)
      dft = kwargs.get('dft', 1)
      
      start_date = kwargs.get('start_date', None)
      self.start_date = pd.to_datetime(start_date)
      
      if 'time_frame' in kwargs:
          time_frame = kwargs.get('time_frame', None)
          self.end_date = self.start_date + pd.to_timedelta(time_frame)
      else:
          end_date = kwargs.get('end_date', None)
          self.end_date = pd.to_datetime(end_date)      
                            
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
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('extracting meteo from {}\n'.format(filename))
        sys.stdout.flush()
        #---------------------------------------------------------------------        

     	try: 
       		f = open(filename)
      	except:
        	print 'no file {}'.format(filename)
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

                
            if (ft1 <= int(stepRange) <= ft2) & (tstamp <= self.end_date):
                
                name,varin,ilon,ilat=getd(gid)    
            
            else:
                grib_release(gid)
                if int(stepRange) > ft2 : break
                continue
                
            lon=ilon[0,:]
            lat=ilat[:,0]
           
        # verbose 
        #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('retrieving {} at {}\n'.format(name, tstamp))
            sys.stdout.flush()
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
    
      self.uvp = met
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
        
        minlon = minlon - .1
        maxlon = maxlon + .1
        minlat = minlat - .1
        maxlat = maxlat + .1
  
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
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('extracting meteo from {}\n'.format(filename))
            sys.stdout.flush()
          #---------------------------------------------------------------------      


            try: 
                f = open(filename)
            except:
                print 'no file {}'.format(filename)
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
        
        #mask non values
        met['msl'] = met.msl.where(met.msl>0)
        met['u10'] = met.u10.where(met.msl>0)
        met['v10'] = met.v10.where(met.msl>0)
                                      
        self.uvp = met
    
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
            print minlon,maxlon,minlat,maxlat
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
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('extracting meteo from {}\n'.format(filename))
            sys.stdout.flush()
        #---------------------------------------------------------------------      

            try: 
                f = open(filename)
            except:
                print 'no file {}'.format(filename)
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
                                              
        self.uvp = met
    
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
              .isel(lon=slice(i0,data.lon.size),lat=slice(j0,j1))
              .sel(time=tslice)
              )
          sh.lon.values = sh.lon.values -360.

          sh1 = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(lon=slice(0,i1),lat=slice(j0,j1))
              .sel(time=tslice)
              )
              
          tot = xr.concat([sh,sh1],dim='lon')
          
      else:            

          tot = (
              data[['prmslmsl','ugrd10m', 'vgrd10m']]
              .isel(lon=slice(i0,i1),lat=slice(j0,j1))
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
        
    
