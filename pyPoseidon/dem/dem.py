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
        if dem == 'gebco08' :
            self.impl = gebco(**kwargs)
        elif dem == 'gebco14' :
            self.impl = gebco(**kwargs)
        elif dem == 'emodnet' :
            self.impl = emodnet(**kwargs)
        else:
            self.impl = erdap(**kwargs)

class gebco08(dem):
    
    def __init__(self,**kwargs):
           
      self.minlon = kwargs.get('minlon', None)
      self.maxlon = kwargs.get('maxlon', None)
      self.minlat = kwargs.get('minlat', None)
      self.maxlat = kwargs.get('maxlat', None)       
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


      minlon = self.minlon #kwargs.get('lon0', {})
      maxlon = self.maxlon #kwargs.get('lon1', {})
      minlat = self.minlat #kwargs.get('lat0', {})
      maxlat = self.maxlat #kwargs.get('lat1', {}) 

      
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
      self.dlats = lats
         
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
       self.ilons = grid_x
       self.ilats = grid_y
    
    def fix(self,**kwargs):  
         
       shpfile = kwargs.get('shoreline', None)
       fix(self,shpfile)

class gebco14(dem):
    
    def __init__(self,**kwargs):
    
      self.minlon = kwargs.get('minlon', None)
      self.maxlon = kwargs.get('maxlon', None)
      self.minlat = kwargs.get('minlat', None)
      self.maxlat = kwargs.get('maxlat', None)       
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


      minlon = self.minlon #kwargs.get('lon0', {})
      maxlon = self.maxlon #kwargs.get('lon1', {})
      minlat = self.minlat #kwargs.get('lat0', {})
      maxlat = self.maxlat #kwargs.get('lat1', {}) 

      
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
      self.dlats = lats
         
      if 'grid_x' in kwargs.keys():
       grid_x = kwargs.get('grid_x', None)
       grid_y = kwargs.get('grid_y', None)
    # interpolate on the given grid
      #flip on lat to make it increasing for RectBivariateSpline
       ilon=lons[0,:]
       ilat=lats[:,0]
    #  sol=scipy.interpolate.RectBivariateSpline(ilon,ilat,topo.T)#,kx=2,ky=2)

       orig = pyresample.geometry.SwathDefinition(lons=self.dlons,lats=self.dlats) # original points
       targ = pyresample.geometry.SwathDefinition(lons=grid_x,lats=grid_y) # target grid

       #mask land
       water = np.ma.masked_array(self.val,self.val>0)
       water.fill_value = 999999
       
       itopo = pyresample.kd_tree.resample_nearest(orig,water,targ,radius_of_influence=50000,fill_value=999999)

    #   itopo=[]
    #   for x,y in zip(grid_x.ravel(),grid_y.ravel()):
    #      itopo.append(sol(x,y).ravel()[0])
    #---------------------------------------------------------------
    #     l=np.abs(ilon-np.float(x)).argmin()
    #     m=np.abs(ilat-np.float(y)).argmin()
    #     xx = ilon[l-1:l+2]
    #     yy = ilat[m-1:m+2]
    #     zz = topo.T[l-1:l+2,m-1:m+2]
    #     fa=scipy.interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
    #     itopo.append(fa(x,y))

    #  itopo=np.array(itopo)
    #  itopo=itopo.reshape(grid_x.shape)

       self.ival = itopo
       self.ilons = grid_x
       self.ilats = grid_y
       
    def fix(self,**kwargs):  
         
       shpfile = kwargs.get('shoreline', None)
       fix(self,shpfile)

class emodnet(dem):

    def __init__(self,**kwargs):
    
      self.minlon = kwargs.get('minlon', None)
      self.maxlon = kwargs.get('maxlon', None)
      self.minlat = kwargs.get('minlat', None)
      self.maxlat = kwargs.get('maxlat', None)       
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

      minlon = self.minlon #kwargs.get('lon0', {})
      maxlon = self.maxlon #kwargs.get('lon1', {})
      minlat = self.minlat #kwargs.get('lat0', {})
      maxlat = self.maxlat #kwargs.get('lat1', {}) 

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
      topo = ncv[n3][j1:j2,i1:i2]
        
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
       
       itopo = pyresample.kd_tree.resample_nearest(orig,self.val,targ,radius_of_influence=50000,fill_value=999999)

       self.ival = itopo
       self.ilons = grid_x
       self.ilats = grid_y
       
       
    def fix(self,**kwargs):  
         
       shpfile = kwargs.get('shoreline', None)
       fix(self,shpfile)

def fix(b,shpfile):
    
    #define coastline
    
    shp = gp.GeoDataFrame.from_file(shpfile)

           
    xp = b.ilons
    yp = b.ilats
    
    
    #put them all in a list
    ls=[]
    for i in range(shp.shape[0]):
        il = shp.loc[i,'geometry']
        try:
            #print len(il)
            for k in range(len(list(il.geoms))): # for MultiLineStrings
               ls.append(list(il.geoms)[k])
        except:
            ls.append(il)
            
    sall = geometry.MultiLineString(ls) # join
    c = ops.linemerge(sall) #merge
    
    #Select the Line Strings that correspond to our grid
    #create a polygon of the grid
    grp=geometry.Polygon([(xp.min(),yp.min()),(xp.min(),yp.max()),(xp.max(),yp.max()),(xp.max(),yp.min())])    
    
    cl=[] #initialize
    #add Polygons if they intersect with the grid
    for i in range(len(c)):
        z = geometry.Polygon(c[i])
        if z.intersects(grp): 
                cl.append(c[i])
    
    cll = geometry.MultiLineString(cl) #join them into a Multiline
    cg = ops.linemerge(cll) #merge parts if possible
    
    gmask = internal(cg,xp, yp)
    
    #check if the grid polygons intersect the shoreline
    gps = []
    igps = []
    for i in range(1,xp.shape[0]):
        for j in range(1,xp.shape[1]):
            p = geometry.Polygon([(xp[i-1,j-1],yp[i-1,j-1]),(xp[i-1,j],yp[i-1,j]),(xp[i,j],yp[i,j]),(xp[i,j-1],yp[i,j-1])])
            if not p.intersects(cg): 
                gps.append(p)
                igps.append([i,j])
                if i-1 == 0 : 
                    igps.append([i-1,j])
                    igps.append([i-1,j-1])
                if j-1 == 0 : 
                    igps.append([i,j-1])
                    igps.append([i-1,j-1])
    
    #create a mask of all cells not intersecting the shoreline
    imask=np.ones(xp.shape, dtype=bool)
    for [i,j] in igps:
        imask[i,j]=0                
    
    # merge gmask and imask
    tmask =  np.logical_and(np.invert(imask),gmask)
        
    #find islands    
    grid = tmask.astype(int).astype(str)
    grid = [list(x) for x in grid]
    
    
    isls = Solution().numIslands(grid)
    #print(isls)
    mgrid = np.array(grid).astype(int)
    nps = []
    for i in range(1, isls + 1):
        nps.append(np.sum(mgrid != i))
    
    val, idx = min((val, idx) for (idx, val) in enumerate(nps))    
    vwater = np.ma.masked_array(mgrid, mgrid != idx + 1)
    vwater[np.invert(vwater.mask)] = 1    
    
    
    fmask = np.logical_and(tmask,np.invert(vwater.mask)) # total mask
    
    ##resample bathymetry
    
    gx=np.ma.masked_array(xp,np.invert(tmask)) #grid to be resampled
    gy=np.ma.masked_array(yp,np.invert(tmask))  
    
    #mask positive bathymetry 
    wet = np.ma.masked_array(b.val,b.val>0)
   # wet.fill_value = 0.
    mx = np.ma.masked_array(b.dlons,b.val>0) 
    my = np.ma.masked_array(b.dlats,b.val>0)
    
    orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=gx,lats=gy) # wet points
    
    # with nearest using only the water values
    b2_near = pyresample.kd_tree.resample_nearest(orig,wet,targ,radius_of_influence=50000,fill_value=999999)
    
    bw = np.ma.masked_array(b2_near,np.invert(fmask))
        
    b.fval = bw
    
    
def internal(cg, xp, yp):
    
    points=zip(xp.flatten(),yp.flatten())
    
    #collect the internal points
    xi=[]
    yi=[]
    for i in range(len(cg)):
        z = geometry.Polygon(cg[i])
        path = mpltPath.Path(zip(z.boundary.xy[0],z.boundary.xy[1]))

        inside = path.contains_points(points)

        if np.sum(inside) > 0:
            X = np.ma.masked_array(xp,mask=np.invert(inside)),
            Y = np.ma.masked_array(yp,mask=np.invert(inside))
            xi.append(X)
            yi.append(Y)    
    
    #merge the masks 
    gmask=np.ones(xi[0][0].shape, dtype=bool)
    for i in range(len(xi)):
        gmask = np.logical_and(gmask,xi[i][0].mask)    
    
    return gmask


class erdap(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
           
      if minlon < 0: minlon = minlon + 360.
      
      if maxlon < 0: maxlon = maxlon + 360.
    
      url = kwargs.get('url', 'http://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm15plus')
            
      data = xr.open_dataset(url)    
      
      i0=np.abs(data.longitude.data-minlon).argmin()
      i1=np.abs(data.longitude.data-maxlon).argmin()

      
      j0=np.abs(data.latitude.data-minlat).argmin()
      j1=np.abs(data.latitude.data-maxlat).argmin()
      
      dem = (
          data.z
          .isel(longitude=slice(i0,i1),latitude=slice(j0,j1))
          )
      
      self.ival = dem
      
      
class gebco(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
           
      if minlon < -180: minlon = minlon + 360.
      
      if maxlon < -180: maxlon = maxlon + 360.
    
      url = kwargs.get('dpath', None)
            
      data = xr.open_dataset(url)    
      
      i0=np.abs(data.lon.data-minlon).argmin()
      i1=np.abs(data.lon.data-maxlon).argmin()

      
      j0=np.abs(data.lat.data-minlat).argmin()
      j1=np.abs(data.lat.data-maxlat).argmin()
      
      dem = (
          data[data.data_vars.keys()[0]]
          .isel(lon=slice(i0,i1),lat=slice(j0,j1))
          )
      
      self.ival = dem
      
class emodnet_(dem):
    
    def __init__(self,**kwargs):
                          
      minlon = kwargs.get('minlon', None)
      maxlon = kwargs.get('maxlon', None)
      minlat = kwargs.get('minlat', None)
      maxlat = kwargs.get('maxlat', None) 
           
    
      url = kwargs.get('dpath', None)
            
      data = xr.open_dataset(url)   
            
      if minlon < data.longitude.min() : 
          print 'longitude out of range [{},{}]'.format(data.longitude[0].values,data.longitude[-1].values)
          sys.exit(1)
      if maxlon > data.longitude.max() : 
          print 'longitude out of range [{},{}]'.format(data.longitude[0].values,data.longitude[-1].values)
          sys.exit(1)
      if minlat < data.latitude.min() : 
          print 'latitude out of range [{},{}]'.format(data.latitude[0].values,data.latitude[-1].values)
          sys.exit(1)
      if maxlat > data.latitude.max() : 
          print 'latitude out of range [{},{}]'.format(data.latitude[0].values,data.latitude[-1].values)
          sys.exit(1)
      
      i0=np.abs(data.longitude.data-minlon).argmin()
      i1=np.abs(data.longitude.data-maxlon).argmin()

      
      j0=np.abs(data.latitude.data-minlat).argmin()
      j1=np.abs(data.latitude.data-maxlat).argmin()
      
      dem = (
          data[data.data_vars.keys()[0]]
          .isel(longitude=slice(i0,i1),latitude=slice(j0,j1))
          )
      
      self.ival = dem

    