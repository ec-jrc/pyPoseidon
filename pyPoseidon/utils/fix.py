"""
Grid adjustment functions

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import geopandas as gp
import matplotlib.path as mpltPath
import shapely 
from pyPoseidon.utils.bfs import *
import pyresample
import pandas as pd
import xarray as xr
import sys
import os
import logging


#logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(levelname)-8s %(asctime)s:%(name)s:%(message)s')

file_handler = logging.FileHandler('dem.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

sformatter = logging.Formatter('%(levelname)-8s %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(sformatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)


def fix(dem,shpfile,nc=50,**kwargs):
    
    #--------------------------------------------------------------------- 
    logger.info('optimize grid\n')
    #--------------------------------------------------------------------- 
    
    ncores = kwargs.get('ncores', 1)
    
    #define coastline
    
    shp = gp.GeoDataFrame.from_file(shpfile)

           
    xp = dem.ilons.values
    yp = dem.ilats.values
    
    
    #put them all in a list
    ls=[]
    for i in range(shp.shape[0]):
        il = shp.loc[i,'geometry']
        try:
            #print(len(il))
            for k in range(len(list(il.geoms))): # for MultiLineStrings
               ls.append(list(il.geoms)[k])
        except:
            ls.append(il)
            
    sall = shapely.geometry.MultiLineString(ls) # join
    c = shapely.ops.linemerge(sall) #merge
    
    #Select the Line Strings that correspond to our grid
    #create a polygon of the grid
    grp=shapely.geometry.Polygon([(xp.min(),yp.min()),(xp.min(),yp.max()),(xp.max(),yp.max()),(xp.max(),yp.min())])    
    
    cl=[] #initialize
    #add Polygons if they intersect with the grid
    for i in range(len(c)):
        z = shapely.geometry.Polygon(c[i])
        if z.intersects(grp): 
                cl.append(c[i])
    
    cll = shapely.geometry.MultiLineString(cl) #join them into a Multiline
    cg = shapely.ops.linemerge(cll) #merge parts if possible
    
    # brake if no coastline
    if cg.geom_type == 'GeometryCollection' :
        wmask = np.zeros(xp.shape, dtype=bool)
        return wmask, cg
    
            
    #--------------------------------------------------------------------- 
    logger.info('..serial..\n')
    #--------------------------------------------------------------------- 

    gps = []
    igps = []
    for i in range(0,xp.shape[0]-1):
        for j in range(0,xp.shape[1]-1):
            p = shapely.geometry.Polygon([(xp[i,j],yp[i,j]),(xp[i+1,j],yp[i+1,j]),(xp[i+1,j+1],yp[i+1,j+1]),(xp[i,j+1],yp[i,j+1])])
            if not p.intersects(cg): 
                gps.append(p)
                igps.append([i+1,j+1])
                if i == 0 : 
                    igps.append([i,j])
                    igps.append([i,j+1])
                if j == 0 : 
                    igps.append([i,j])
                    igps.append([i+1,j])

    ps = igps                
    #------------------------------------------------------------------------------
        
    #create a mask of all cells not intersecting the shoreline
    imask=np.ones(xp.shape, dtype=bool)
    for [i,j] in ps:
        imask[i,j]=0                
    
    #Find points inside coastline (land)
    gmask = internal(cg,xp,yp)
 
    # merge gmask and imask
    tmask =  np.logical_and(np.invert(imask),gmask)
    tmask = np.invert(tmask)         
    
    # break if there are no internal points
    if gmask.sum() == 0 : 
        return imask, cg
        
    #--------------------------------------------------------------------- 
    logger.info('eliminate isolated wet regions\n')
    #--------------------------------------------------------------------- 
                
    #find islands    
    p = rmt(tmask,xp,nc) 
    
    # repeat until convergence
    k=0
    while True:
        p1 = rmt(p,xp,nc)
        k+=1
        #--------------------------------------------------------------------- 
        logger.info('... adjusting ...\n')
        #---------------------------------------------------------------------  
        if np.array_equal(p1,p) : break
        p = p1
    
    wmask = p1   
    
    #--------------------------------------------------------------------- 
    logger.info('done \n')
    #--------------------------------------------------------------------- 
    
    return wmask, cg
        
    
def bmatch(dem,wmask,**kwargs):
       
    #--------------------------------------------------------------------- 
    logger.info('resample bathymetry\n')
    #--------------------------------------------------------------------- 
    
    ncores = kwargs.get('ncores', 1)
       
    xp = dem.ilons.values
    yp = dem.ilats.values
          
    ##resample bathymetry
    xw=np.ma.masked_array(xp,wmask) #wet points 
    yw=np.ma.masked_array(yp,wmask) 
    
    # fill the nan, if present, we values in order to compute values there if needed.
    dem.elevation.data[np.isnan(dem.elevation)]=9999. 
    
    #mask positive bathymetry 
    wet = np.ma.masked_array(dem.elevation,dem.elevation>0)
    x, y = np.meshgrid(dem.longitude,dem.latitude)
   # wet.fill_value = 0.
    mx = np.ma.masked_array(x,dem.elevation.values>0) 
    my = np.ma.masked_array(y,dem.elevation.values>0)
    
    orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=xw,lats=yw) # wet points
    
    # with nearest using only the water values
    
    grid_con = pyresample.image.ImageContainerNearest(wet, orig, radius_of_influence=50000,fill_value=np.nan)#, nprocs=ncores)

    area_con = grid_con.resample(targ)

    result = area_con.image_data
                   
    bw = np.ma.masked_array(result,wmask)
        
    if np.isnan(bw).all() == True: # use the ivals if all is Nan, e.g. only land
        bw = dem.ival.values
        
    fval = xr.Dataset({'fval': (['ilat', 'ilon'],  bw)},   
                                  coords={'ilon': ('ilon', xp[0,:]),   
                                          'ilat': ('ilat', yp[:,0])})
    
    dem = xr.merge([dem,fval])
    
    #--------------------------------------------------------------------- 
    logger.info('done \n')
    #--------------------------------------------------------------------- 
    
    return dem
    
def internal(cg, xp, yp):
    
    points=list(zip(xp.flatten(),yp.flatten()))
    
    #collect the internal points
    xi=[]
    yi=[]
    
    if cg.geom_type.startswith('Multi'): # many LineStrings
    
        for i in range(len(cg)):
            z = shapely.geometry.Polygon(cg[i])
            path = mpltPath.Path(list(zip(z.boundary.xy[0],z.boundary.xy[1])))

            inside = path.contains_points(points)

            if np.sum(inside) > 0:
                X = np.ma.masked_array(xp,mask=np.invert(inside)),
                Y = np.ma.masked_array(yp,mask=np.invert(inside))
                xi.append(X)
                yi.append(Y) 
                
    else: #Single LineString
    
        z = shapely.geometry.Polygon(cg)
        path = mpltPath.Path(list(zip(z.boundary.xy[0],z.boundary.xy[1])))

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
    for key, val in mask.items():
        fmask = np.logical_xor(val.mask,fmask)

    if ne % 2 == 0:    
        return np.invert(fmask) # this is wet
    else:
        return fmask



def loop(xp,yp,cg):
#check if the grid polygons intersect the shoreline
    gps = []
    igps = []
    for i in range(xp.shape[0]-1):
        for j in range(xp.shape[1]-1):
            p = shapely.geometry.Polygon([(xp[i,j],yp[i,j]),(xp[i+1,j],yp[i+1,j]),(xp[i+1,j+1],yp[i+1,j+1]),(xp[i,j+1],yp[i,j+1])])
            if not p.intersects(cg): 
                gps.append(p)
                igps.append([i+1,j+1])
    
    return igps