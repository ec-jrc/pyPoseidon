"""
Grid adjustment functions

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import geopandas as gp
import shapely 
import pygeos
from pyPoseidon.utils.bfs import *
import pyresample
import pandas as pd
import xarray as xr
import sys
import os
import logging


#logging setup
import logging
logger = logging.getLogger('pyPoseidon')

def fix(dem,coastline,**kwargs):
    
    #--------------------------------------------------------------------- 
    logger.info('adjust dem\n')
    #--------------------------------------------------------------------- 

    #define coastline    
    try:
        shp = gp.GeoDataFrame.from_file(coastline)
    except:
        shp = gp.GeoDataFrame(coastline)

    
    if 'ival' in dem.data_vars:       
        xp = dem.ilons.values
        yp = dem.ilats.values
    else:
        xp = dem.longitude.values
        yp = dem.latitude.values
        
    minlon = xp.min()
    maxlon = xp.max()
    minlat = yp.min()
    maxlat = yp.max()

    
    if xp.mean() < 0 and xp.min() < -180. :
        flag = -1
    elif xp.mean() > 0 and xp.max() > 180. :
        flag = 1
    else:
        flag = 0
    
    
    if flag == 1 :
        block1 = shp.cx[minlon:180,minlat:maxlat].copy()
        block2 = shp.cx[-180:(maxlon-360.),minlat:maxlat].copy()
    
        for idx, poly in block2.iterrows():
            block2.loc[idx,'geometry'] = shapely.ops.transform(lambda x,y,z=None: (x + 360.,y), poly.geometry)
    
        block = block1.append(block2)
    
    elif flag == -1 :
    
        block1 = shp.cx[minlon + 360 : 180,minlat:maxlat].copy()
        block2 = shp.cx[-180:maxlon,minlat:maxlat].copy()
    
        for idx, poly in block1.iterrows():
            block1.loc[idx,'geometry'] = shapely.ops.transform(lambda x,y,z=None: (x - 360.,y), poly.geometry)
    
        block = block1.append(block2)

    else:
        block = shp.cx[minlon:maxlon,minlat:maxlat]
        
    try:
        block = gp.GeoDataFrame(geometry = list(block.unary_union))
    except:
        pass

    #--------------------------------------------------------------------- 
    logger.debug('compute water and land\n')
    #--------------------------------------------------------------------- 

    
    #create a polygon of the lat/lon window
    grp=shapely.geometry.Polygon([(minlon,minlat),(minlon,maxlat),(maxlon,maxlat),(maxlon,minlat)])

    grp = grp.buffer(.5) # buffer it to get also the boundary points
    
    g = block.unary_union.symmetric_difference(grp) # get the diff
    
    try:
        t = gp.GeoDataFrame({'geometry':g})
    except:
        t = gp.GeoDataFrame({'geometry':[g]})

    t['length']=t['geometry'][:].length # optional
    
    t = t.sort_values(by='length', ascending=0) #use the length to list them
    t = t.reset_index(drop=True)
    
    t['in'] = gp.GeoDataFrame(geometry=[grp.buffer(.001)] * t.shape[0]).contains(t) # find the largest of boundaries
    
    try:
        idx = np.where(t['in']==True)[0][0] # first(largest) boundary within lat/lon
        b = t.iloc[idx].geometry #get the largest 
    except:
        b = shapely.geometry.GeometryCollection()

    #define wet/dry
    water = b
    land = grp - b
    
    if (not land) | (not water):
        
        #--------------------------------------------------------------------- 
        logger.debug('only water/land present...\n')
        #--------------------------------------------------------------------- 
        
        if 'ival' in dem.data_vars:
                
            dem = dem.assign(fval=dem.ival)
        
        else:
        
            dem = dem.assign(adjusted=dem.elevation)
            
        return dem
    
    if 'ival' in dem.data_vars:
        df = pd.DataFrame({'longitude':dem.ilons.values.flatten(),'latitude':dem.ilats.values.flatten(),'elevation':dem.ival.values.flatten()})
    else:
        df = dem.elevation.to_dataframe().reset_index()
    
    #--------------------------------------------------------------------- 
    logger.debug('invoke pygeos\n')
    #--------------------------------------------------------------------- 
    
    spoints_ = pygeos.points(list(df.loc[:,['longitude','latitude']].values)) # create pygeos objects for the points
    
    # Add land boundaries to a pygeos object
    try:
        lbs = []
        for l in range(len(land.boundary)):
            z = pygeos.linearrings(land.boundary[l].coords[:]) 
            lbs.append(z)        
    except:
        lbs = pygeos.linearrings(land.boundary.coords[:])
        
            
    bp = pygeos.polygons(lbs)
    
    #--------------------------------------------------------------------- 
    logger.debug('find wet and dry masks\n')
    #--------------------------------------------------------------------- 
    
    #find the points on land
    
    tree = pygeos.STRtree(spoints_)
    
    try:
        wl=[]
        for l in range(len(land.boundary)):
            wl.append(tree.query(bp[l], predicate='contains').tolist())
    except:
        wl = tree.query(bp, predicate='contains').tolist()
    
    ns = [j for i in wl for j in i]
    
    lmask=np.zeros(spoints_.shape, dtype=bool)
    lmask[ns] = True
       
    wmask = ~lmask # invert for wet mask
    
    #--------------------------------------------------------------------- 
    logger.debug('fix wet points\n')
    #--------------------------------------------------------------------- 
    
    
    #Now see if the wet points have indeed negative values
    
    pw_mask = df.loc[wmask,'elevation'] > 0 
    
    if pw_mask.sum() > 0:
        
        pw = df.loc[wmask][pw_mask] # problematic points: bathymetry > 0 in wet area
    
        #Resample to fix that ...
        xw = pw.longitude.values
        yw = pw.latitude.values
        
        #Define points with positive bathymetry        
        x, y = np.meshgrid(dem.longitude,dem.latitude)#!!!!!!!!
        
        if flag == 1:
            xw = xw - 180.
            x = x - 180.
        elif flag == -1:
            xw = xw + 180.
            x = x + 180.
        
        # wet.fill_value = 0.
        mx = np.ma.masked_array(x,dem.elevation.values>0) 
        my = np.ma.masked_array(y,dem.elevation.values>0)
    
        # fill the nan, if present, with values in order to compute values there if needed.
        dem.elevation.data[np.isnan(dem.elevation.values)]=9999. 

        #mask positive bathymetry 
        wet_dem = np.ma.masked_array(dem.elevation,dem.elevation.values > 0 )
    
        orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original bathymetry points
        targ = pyresample.geometry.SwathDefinition(lons=xw,lats=yw) # wet points
        
        bw = pyresample.kd_tree.resample_nearest(orig,wet_dem,targ,radius_of_influence=50000,fill_value=np.nan)
        
        df.loc[pw.index,'elevation'] = bw # replace in original dataset
    
    
    #--------------------------------------------------------------------- 
    logger.debug('fix dry points\n')
    #--------------------------------------------------------------------- 
    
    #.. the same for dry points
    
    pl_mask = df.loc[lmask,'elevation'] < 0 
    
    if pl_mask.sum() > 0:
        
        pl = df.loc[lmask][pl_mask] # problematic points: bathymetry <0 in dry area
    
        ## Resample to fix that 
        xl = pl.longitude.values
        yl = pl.latitude.values
        
        x, y = np.meshgrid(dem.longitude,dem.latitude)
        
        if flag == 1:
            xl = xl - 180.
            x = x - 180.
        elif flag == -1:
            xl = xl + 180.
            x = x + 180.
        
        # wet.fill_value = 0.
        dx = np.ma.masked_array(x,dem.elevation.values<0) 
        dy = np.ma.masked_array(y,dem.elevation.values<0)
        
        # fill the nan, if present, with values in order to compute values there if needed.
        dem.elevation.data[np.isnan(dem.elevation.values)]=9999. 

        #mask positive bathymetry 
        dry_dem = np.ma.masked_array(dem.elevation,dem.elevation.values < 0 )
    
        orig = pyresample.geometry.SwathDefinition(lons=dx,lats=dy) # original bathymetry points
        targ = pyresample.geometry.SwathDefinition(lons=xl,lats=yl) # wet points
        
        bd = pyresample.kd_tree.resample_nearest(orig,dry_dem,targ,radius_of_influence=50000,fill_value=np.nan)
        
        df.loc[pl.index,'elevation'] = bd  # replace in original dataset
        
    #--------------------------------------------------------------------- 
    logger.debug('assemble dataset \n')
    #--------------------------------------------------------------------- 
    
    
    #reassemble dataset
    
    if 'ival' in dem.data_vars:
        
        if len(dem.ival.shape) == 1 :
            new_dem = df.elevation.to_xarray()
            new_dem = xr.merge([new_dem])
            new_dem = new_dem.rename({'elevation':'fval'})
            new_dem.fval.attrs = {'coastline':'based on coastline'}
            new_dem = new_dem.rename({'index':'k'}).drop_vars('k')
        else:
            new_dem = df.set_index(['latitude','longitude']).to_xarray().rename({'longitude':'l','latitude':'k','elevation':'fval'}).drop_vars(['k','l'])

        
    else:
        
        df_new = df.set_index(['latitude','longitude'])
        new_dem = df_new.to_xarray() 
        new_dem = new_dem.rename({'elevation':'adjusted'})
        new_dem.attrs = {'coastline':'based on coastline'}
    
    cdem = xr.merge([dem,new_dem])
    
    return cdem

