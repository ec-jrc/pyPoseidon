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

try:
    import pp
except:
    pass

def fix(dem,shpfile,**kwargs):
    
    #--------------------------------------------------------------------- 
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('optimize grid\n')
    sys.stdout.flush()
    #--------------------------------------------------------------------- 
    
    nc = kwargs.get('nc', 10)
    ncores = kwargs.get('ncores', 1)
    
    #define coastline
    
    shp = gp.GeoDataFrame.from_file(shpfile)

           
    xp = dem.Dataset.ilons.values
    yp = dem.Dataset.ilats.values
    
    
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
    
    
    try:
        #------------------------------------------------------------------------------
        #check if the grid polygons intersect the shoreline USING PP
        try:
            PYTHONPATH =  os.environ['PYTHONPATH'] #SAVE PYTHONPATH in order to reset it afterwards
        except:
            pass
    
        job_server = pp.Server() 
       
        if ncores > job_server.get_ncpus(): ncores = job_server.get_ncpus() # make sure we don't overclock
    
        job_server.set_ncpus(ncores)
    
        #split the array
    
        n = xp.shape[0]
        l=range(n)
        k = ncores
    
        b = [l[i * (n // k) + min(i, n % k):(i+1) * (n // k) + min(i+1, n % k)] for i in range(k)]
    
        jobs=[]
        for l in range(ncores):
            extra=b[l][-1]+1 # add one to consider the jump between parts
            part = b[l]+[extra]
            if part[-1] >= xp.shape[0] : part = part[:-1] #control the last element
            xl=xp[part,:]
            yl=yp[part,:]   

            jobs.append(job_server.submit(loop,(xl,yl,cg,),modules=('shapely',)))

        ps=jobs[0]()
        for l in range(1,ncores):    
            ps = np.vstack([ps, [[i+b[l][0],j] for [i,j] in jobs[l]()]]) # fix global index by adding the first index of part    
    
        #Add extra points from boundary    
        bgps=[]
        for [i,j] in ps:
            if i == 1 : 
                bgps.append([i-1,j-1])
                bgps.append([i-1,j])
            if j == 1 : 
                bgps.append([i-1,j-1])
                bgps.append([i,j-1])
            
        ps = np.vstack([ps, bgps]) # final stack
    
        job_server.destroy()
    
        try:
            if PYTHONPATH :
                os.environ['PYTHONPATH'] = PYTHONPATH  #reset PYTHONPATH 
            else:    
                try :
                    del os.environ['PYTHONPATH']
                except:
                    pass 
        except:
            pass
            
    except :
        
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('..serial..\n')
        sys.stdout.flush()
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
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('eliminate isolated wet regions\n')
    sys.stdout.flush()
    #--------------------------------------------------------------------- 
                
    #find islands    
    p = rmt(tmask,xp,nc) 
    
    # repeat until convergence
    k=0
    while True:
        p1 = rmt(p,xp,nc)
        k+=1
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('... adjusting ...\n')
        sys.stdout.flush()
        #---------------------------------------------------------------------  
        if np.array_equal(p1,p) : break
        p = p1
    
    wmask = p1   
    
    #--------------------------------------------------------------------- 
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('done \n')
    sys.stdout.flush()
    #--------------------------------------------------------------------- 
    
    return wmask, cg
        
    
def bmatch(dem,wmask,**kwargs):
       
    #--------------------------------------------------------------------- 
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('resample bathymetry\n')
    sys.stdout.flush()
    #--------------------------------------------------------------------- 
    
    ncores = kwargs.get('ncores', 1)
       
    xp = dem.Dataset.ilons.values
    yp = dem.Dataset.ilats.values
          
    ##resample bathymetry
    xw=np.ma.masked_array(xp,wmask) #wet points 
    yw=np.ma.masked_array(yp,wmask) 
    
    # fill the nan, if present, we values in order to compute values there if needed.
    dem.Dataset.val.data[np.isnan(dem.Dataset.val)]=9999. 
    
    #mask positive bathymetry 
    wet = np.ma.masked_array(dem.Dataset.val,dem.Dataset.val>0)
   # wet.fill_value = 0.
    mx = np.ma.masked_array(dem.Dataset.dlons,dem.Dataset.val.values>0) 
    my = np.ma.masked_array(dem.Dataset.dlats,dem.Dataset.val.values>0)
    
    orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original bathymetry points
    targ = pyresample.geometry.SwathDefinition(lons=xw,lats=yw) # wet points
    
    # with nearest using only the water values
    
    grid_con = pyresample.image.ImageContainerNearest(wet, orig, radius_of_influence=50000,fill_value=np.nan)#, nprocs=ncores)

    area_con = grid_con.resample(targ)

    result = area_con.image_data
                   
    bw = np.ma.masked_array(result,wmask)
        
    if np.isnan(bw).all() == True: # use the ivals if all is Nan, e.g. only land
        bw = dem.Dataset.ival.values
        
    fval = xr.Dataset({'fval': (['ilat', 'ilon'],  bw)},   
                                  coords={'ilon': ('ilon', xp[0,:]),   
                                          'ilat': ('ilat', yp[:,0])})
    
    dem.Dataset = xr.merge([dem.Dataset,fval])
    
    #--------------------------------------------------------------------- 
    sys.stdout.flush()
    sys.stdout.write('\n')
    sys.stdout.write('done \n')
    sys.stdout.flush()
    #--------------------------------------------------------------------- 
    
    
    
def internal(cg, xp, yp):
    
    points=zip(xp.flatten(),yp.flatten())
    
    #collect the internal points
    xi=[]
    yi=[]
    
    if cg.geom_type.startswith('Multi'): # many LineStrings
    
        for i in range(len(cg)):
            z = shapely.geometry.Polygon(cg[i])
            path = mpltPath.Path(zip(z.boundary.xy[0],z.boundary.xy[1]))

            inside = path.contains_points(points)

            if np.sum(inside) > 0:
                X = np.ma.masked_array(xp,mask=np.invert(inside)),
                Y = np.ma.masked_array(yp,mask=np.invert(inside))
                xi.append(X)
                yi.append(Y) 
                
    else: #Single LineString
    
        z = shapely.geometry.Polygon(cg)
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