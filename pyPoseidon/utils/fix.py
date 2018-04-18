import numpy as np
import geopandas as gp
import matplotlib.path as mpltPath
from shapely import geometry, ops
from mpl_toolkits.basemap import Basemap
from pyPoseidon.utils.bfs import *
import pyresample
import pandas as pd
import xarray as xr

def fix(b,shpfile, nc=10):
    
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
    
    
    #check if the grid polygons intersect the shoreline
    gps = []
    igps = []
    for i in range(0,xp.shape[0]-1):
        for j in range(0,xp.shape[1]-1):
            p = geometry.Polygon([(xp[i,j],yp[i,j]),(xp[i+1,j],yp[i+1,j]),(xp[i+1,j+1],yp[i+1,j+1]),(xp[i,j+1],yp[i,j+1])])
            if not p.intersects(cg): 
                gps.append(p)
                igps.append([i+1,j+1])
                if i == 0 : 
                    igps.append([i,j])
                    igps.append([i,j+1])
                if j == 0 : 
                    igps.append([i,j])
                    igps.append([i+1,j])
    
    #create a mask of all cells not intersecting the shoreline
    imask=np.ones(xp.shape, dtype=bool)
    for [i,j] in igps:
        imask[i,j]=0                
    
    #Find points inside coastline (land)
    gmask = internal(cg,xp,yp)
    
    
    # merge gmask and imask
    tmask =  np.logical_and(np.invert(imask),gmask)
    tmask = np.invert(tmask) 
        
    #find islands    
    p = rmt(tmask,xp,nc) 
    
    # repeat until convergence
    k=0
    while True:
        p1 = rmt(p,xp,nc)
        k+=1
        print k
        if np.array_equal(p1,p) : break
        p = p1
    
    wmask = p1
    
       
    ##resample bathymetry
    xw=np.ma.masked_array(xp,wmask) #wet points 
    yw=np.ma.masked_array(yp,wmask) 
    
    #mask positive bathymetry 
    wet = np.ma.masked_array(b.val,b.val>0)
   # wet.fill_value = 0.
    mx = np.ma.masked_array(b.dlons,b.val>0) 
    my = np.ma.masked_array(b.dlats,b.val>0)
    
    orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=xw,lats=yw) # wet points
    
    # with nearest using only the water values
    b2_near = pyresample.kd_tree.resample_nearest(orig,wet,targ,radius_of_influence=50000,fill_value=999999)
    
    bw = np.ma.masked_array(b2_near,wmask)
        
    fval = xr.Dataset({'fval': (['latitude', 'longitude'],  bw)},   
                                  coords={'longitude': ('longitude', xp[0,:]),   
                                          'latitude': ('latitude', yp[:,0])})
    
    b.dem = xr.merge([b.dem,fval])
    
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

def rmt(imask,xp,nc):
    grid = np.invert(imask).astype(int).astype(str)
    grid = [list(x) for x in grid]
    isls = Solution().numIslands(grid)

    mgrid = np.array(grid).astype(int)
    nps = []
    for i in range(1, isls + 1):
        nps.append(np.sum(mgrid != i))
    
    grs = pd.DataFrame([(val, idx) for (idx, val) in enumerate(nps)],columns=['val','idx'])

    grs = grs.sort_values('val').reset_index(drop=True)
    
    grs['np']=mgrid.size-grs.val # count number of points in islands
    
    ne = grs[grs.np>nc].index.max() + 1
        
    mask={}
    for i in range(ne):
        b = mgrid.copy()
        vm = np.ma.masked_array(b, b != grs.loc[i,'idx']+1)
        vm[np.invert(vm.mask)] = 1    
        mask.update({str(i):vm})

    fmask=np.zeros(xp.shape, dtype=bool) # initiate

#merge masks
    for key, val in mask.iteritems():
        fmask = np.logical_xor(val.mask,fmask)

    if ne % 2 == 0:    
        return np.invert(fmask) # this is wet
    else:
        return fmask

