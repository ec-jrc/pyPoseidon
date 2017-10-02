import numpy as np
import netCDF4
import scipy.interpolate
import pyresample

class dem:
    impl=None
    def __init__(self, **kwargs):
        dem = kwargs.get('dem', None)
        if dem == 'gebco08' :
            self.impl = gebco08(**kwargs)
        elif dem == 'gebco14' :
            self.impl = gebco14(**kwargs)
        elif dem == 'emodnet' :
            self.impl = emodnet(**kwargs)


class gebco08(dem):
    
    def __init__(self,**kwargs):
           
      self.lon0 = kwargs.get('lon0', None)
      self.lon1 = kwargs.get('lon1', None)
      self.lat0 = kwargs.get('lat0', None)
      self.lat1 = kwargs.get('lat1', None)       
      self.properties = kwargs.get('properties', {})
      
      filename = kwargs.get('dpath', None) 
                
    # open NetCDF data in 
      nc = netCDF4.Dataset(filename)
      ncv = nc.variables
    #print ncv.keys()
      nv=ncv.keys()
     
      n1,n2,n3=nv[0],nv[1],nv[2] # lon,lat,elevation

      lon = ncv[n1][:]
      lat = ncv[n2][:]


      minlon = self.lon0 #kwargs.get('lon0', {})
      maxlon = self.lon1 #kwargs.get('lon1', {})
      minlat = self.lat0 #kwargs.get('lat0', {})
      maxlat = self.lat1 #kwargs.get('lat1', {}) 

      
      if minlon > 175:

        lon=lon+180.

        i1=np.abs(lon-minlon).argmin() 
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3

        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

        zlon=lon.shape[0]

        topo = ncv[n3][j1:j2,zlon/2+i1:]
        topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

      elif minlon < -175:

        lon=lon-180.

        i1=np.abs(lon-minlon).argmin()
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3


        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

        zlon=lon.shape[0]

        topo = ncv[n3][j1:j2,zlon/2+i1:]
        topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

      else:

        i1=np.abs(lon-minlon).argmin()
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3

        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
        topo = ncv[n3][j1:j2,i1:i2]
        
      self.val = topo
      self.dlons = lons
      self.dlat = lats
         
      if 'grid_x' in kwargs.keys():
       grid_x = kwargs.get('grid_x', None)
       grid_y = kwargs.get('grid_y', None)
    # interpolate on the given grid
      #flip on lat to make it increasing for RectBivariateSpline
       ilon=lons[0,:]
       ilat=lats[:,0]
       sol=scipy.interpolate.RectBivariateSpline(ilon,ilat,topo.T,kx=2,ky=2)

       itopo=[]
       for x,y in zip(grid_x.ravel(),grid_y.ravel()):
          itopo.append(sol(x,y).ravel()[0])
    #---------------------------------------------------------------
    #     l=np.abs(ilon-np.float(x)).argmin()
    #     m=np.abs(ilat-np.float(y)).argmin()
    #     xx = ilon[l-1:l+2]
    #     yy = ilat[m-1:m+2]
    #     zz = topo.T[l-1:l+2,m-1:m+2]
    #     fa=scipy.interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
    #     itopo.append(fa(x,y))

       itopo=np.array(itopo)
       itopo=itopo.reshape(grid_x.shape)

       self.ival = itopo


class gebco14(dem):
    
    def __init__(self,**kwargs):
    
      self.lon0 = kwargs.get('lon0', None)
      self.lon1 = kwargs.get('lon1', None)
      self.lat0 = kwargs.get('lat0', None)
      self.lat1 = kwargs.get('lat1', None)       
      self.properties = kwargs.get('properties', {})    
        
      filename = kwargs.get('dpath', None)      
    # open NetCDF data in 
      nc = netCDF4.Dataset(filename)
      ncv = nc.variables
    #print ncv.keys()
      nv=ncv.keys()
     
      n3,n2,n1=nv[0],nv[1],nv[2]

      lon = ncv[n1][:]
      lat = ncv[n2][:]


      minlon = self.lon0 #kwargs.get('lon0', {})
      maxlon = self.lon1 #kwargs.get('lon1', {})
      minlat = self.lat0 #kwargs.get('lat0', {})
      maxlat = self.lat1 #kwargs.get('lat1', {}) 

      
      if minlon > 175:

        lon=lon+180.

        i1=np.abs(lon-minlon).argmin() 
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3

        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

        zlon=lon.shape[0]

        topo = ncv[n3][j1:j2,zlon/2+i1:]
        topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

      elif minlon < -175:

        lon=lon-180.

        i1=np.abs(lon-minlon).argmin()
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3


        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])

        zlon=lon.shape[0]

        topo = ncv[n3][j1:j2,zlon/2+i1:]
        topo = np.hstack([topo,ncv[n3][j1:j2,:i2-zlon/2]])

      else:

        i1=np.abs(lon-minlon).argmin()
        if lon[i1] > minlon: i1=i1-3
        i2=np.abs(lon-maxlon).argmin()
        if lon[i2] < maxlon: i2=i2+3

        j1=np.abs(lat-minlat).argmin()
        if lat[j1] > minlat: j1=j1-3
        j2=np.abs(lat-maxlat).argmin()
        if lat[j2] < maxlat: j2=j2+3

        lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
        topo = ncv[n3][j1:j2,i1:i2]
        
      self.val = topo
      self.dlons = lons 
      self.dlat = lats
         
      if 'grid_x' in kwargs.keys():
       grid_x = kwargs.get('grid_x', None)
       grid_y = kwargs.get('grid_y', None)
    # interpolate on the given grid
      #flip on lat to make it increasing for RectBivariateSpline
       ilon=lons[0,:]
       ilat=lats[:,0]
       sol=scipy.interpolate.RectBivariateSpline(ilon,ilat,topo.T)#,kx=2,ky=2)

       itopo=[]
       for x,y in zip(grid_x.ravel(),grid_y.ravel()):
          itopo.append(sol(x,y).ravel()[0])
    #---------------------------------------------------------------
    #     l=np.abs(ilon-np.float(x)).argmin()
    #     m=np.abs(ilat-np.float(y)).argmin()
    #     xx = ilon[l-1:l+2]
    #     yy = ilat[m-1:m+2]
    #     zz = topo.T[l-1:l+2,m-1:m+2]
    #     fa=scipy.interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
    #     itopo.append(fa(x,y))

       itopo=np.array(itopo)
       itopo=itopo.reshape(grid_x.shape)

       self.ival = itopo


class emodnet(dem):

    def __init__(self,**kwargs):
    
      self.lon0 = kwargs.get('lon0', None)
      self.lon1 = kwargs.get('lon1', None)
      self.lat0 = kwargs.get('lat0', None)
      self.lat1 = kwargs.get('lat1', None)       
      self.properties = kwargs.get('properties', {})    
        
      filename = kwargs.get('dpath', None)      
    # open NetCDF data in 
      nc = netCDF4.Dataset(filename)
      ncv = nc.variables
    #print ncv.keys()
      nv=ncv.keys()
     
      n1,n2,n3=nv[0],nv[1],nv[2]

      lon = ncv[n1][:]
      lat = ncv[n2][:]

      minlon = self.lon0 #kwargs.get('lon0', {})
      maxlon = self.lon1 #kwargs.get('lon1', {})
      minlat = self.lat0 #kwargs.get('lat0', {})
      maxlat = self.lat1 #kwargs.get('lat1', {}) 

      if (minlon < lon.min()) or (maxlon > lon.max()): print 'Lon must be within {} and {}'.format(lon.min(),lon.max())
      if (minlat < lat.min()) or (maxlat > lat.max()): print 'Lat must be within {} and {}'.format(lat.min(),lat.max())

      i1=np.abs(lon-minlon).argmin()
      if lon[i1] > minlon: i1=i1-3
      i2=np.abs(lon-maxlon).argmin()
      if lon[i2] < maxlon: i2=i2+3

      j1=np.abs(lat-minlat).argmin()
      if lat[j1] > minlat: j1=j1-3
      j2=np.abs(lat-maxlat).argmin()
      if lat[j2] < maxlat: j2=j2+3

      lons, lats = np.meshgrid(lon[i1:i2],lat[j1:j2])
      topo = ncv[n3][j1:j2,i1:i2]
        
      self.val = topo
      self.dlons = lons 
      self.dlat = lats
         
      if 'grid_x' in kwargs.keys():
       grid_x = kwargs.get('grid_x', None)
       grid_y = kwargs.get('grid_y', None)
    # resample on the given grid
       ilon=lons[0,:]
       ilat=lats[:,0]
              
       orig = pyresample.geometry.SwathDefinition(lons=self.dlons,lats=self.dlat) # original points
       targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid
       
       itopo = pyresample.kd_tree.resample_nearest(orig,self.val,targ,radius_of_influence=50000,fill_value=999999)
       
       self.ival = itopo


