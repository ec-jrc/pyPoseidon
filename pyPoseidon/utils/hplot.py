"""
Visualization module based on holoviews

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 


import numpy as np
import holoviews as hv
import geoviews as gv
from holoviews import opts
from holoviews.operation.datashader import datashade, rasterize
from datashader.colors import viridis
import geoviews.feature as gf
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import pandas as pd
import panel as pn
from hvplot import xarray
from pyPoseidon.utils.quads2tr import quads_to_tris
gv.extension('bokeh')

import sys
import os
         

@xr.register_dataset_accessor('hplot')
#@xr.register_dataarray_accessor('pplot')

class hplot(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj    
        
     
    def contourf(self, var='depth', it=None, **kwargs):
                
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        try:
            t = kwargs.get('t',self._obj.time.values)
        except:
            pass
        tes = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:4])
        
        # sort out quads
        mask = np.isnan(tes)[:,3]      
        tr3 = tes[mask][:,:3]
        tr3_ = quads_to_tris(tes[~mask])
        tri3 = np.append(tr3,tr3_,axis=0).astype(int)
               
        z = kwargs.get('z',self._obj[var].values[it,:].flatten())
        
        nodes = pd.DataFrame({'longitude':x,'latitude':y, '{}'.format(var):z})
        elems  = pd.DataFrame(tri3, columns=['a', 'b', 'c'])
                    
        width = kwargs.get('width', 800)
        height = kwargs.get('height', 600)        
        opts.defaults(opts.WMTS(width=width, height=height))

        tiles = gv.WMTS('https://maps.wikimedia.org/osm-intl/{Z}/{X}/{Y}@2x.png')

        points = gv.operation.project_points(gv.Points(nodes, vdims=['{}'.format(var)]))

        trimesh = gv.TriMesh((elems, points))
        return tiles * rasterize(trimesh, aggregator='mean').opts(colorbar=True, cmap='Viridis', padding=.1, tools=['hover'])
            
    
    def meteo(self,**kwargs):
        
        m = self._obj.set_coords(['lon','lat'])
        m['lon'] = m['lon'].isel(time=1, drop=True)
        m['lat'] = m['lat'].isel(time=1, drop=True)
        
        return m.prmsl.hvplot.quadmesh(x='lon',y='lat',geo=True,  widget_location='bottom', datashade=True, tiles=True)
    
    def quiver(self,**kwargs):
                
        return   
        
    def grid(self, **kwargs):
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        tes = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values)
        
        width = kwargs.get('width', 800)
        height = kwargs.get('height', 600)        
        opts.defaults(opts.WMTS(width=width, height=height))
        
        tiles = gv.WMTS('https://maps.wikimedia.org/osm-intl/{Z}/{X}/{Y}@2x.png')
        
        nodes = pd.DataFrame({'longitude':x,'latitude':y})        
        elems  = pd.DataFrame(tes, columns=['a', 'b', 'c', 'd'])

        points = gv.operation.project_points(gv.Points(nodes))
        
        if elems.d.is_unique :
            
            trimesh=gv.TriMesh((elems, points)).edgepaths
            
            return tiles * datashade(trimesh, precompute=True, cmap=['black'])        
            
            
        else: # there are quads
        
            quads = elems.loc[~elems.d.isna()].copy()
            quads = quads.reset_index(drop=True)
            ap = nodes.loc[quads.a,['longitude','latitude']]
            bp = nodes.loc[quads.b,['longitude','latitude']]
            cp = nodes.loc[quads.c,['longitude','latitude']]
            dp = nodes.loc[quads.d,['longitude','latitude']]
            quads['ap'] = ap.values.tolist()
            quads['bp'] = bp.values.tolist()
            quads['cp'] = cp.values.tolist()
            quads['dp'] = dp.values.tolist()

            q = gv.Polygons([quads.loc[i,['ap','bp','cp','dp']].tolist() for i in range(quads.shape[0])]).options(fill_alpha=0, line_color='black')
            
            triangles = elems.loc[elems.d.isna()]
            
            trimesh=gv.TriMesh((triangles, points)).edgepaths
            
            g1 = tiles * datashade(trimesh, precompute=True, cmap=['black']) 
            
            
            return g1 * q        
        
    
    def qframes(self,**kwargs):
        
        
        return 
                
 
    def frames(self,var='depth',**kwargs):

        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        try:
            t = kwargs.get('t',self._obj.time.values)
        except:
            pass
        tes = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:4])
        
        # sort out quads
        mask = np.isnan(tes)[:,3]      
        tr3 = tes[mask][:,:3]
        tr3_ = quads_to_tris(tes[~mask])
        tri3 = np.append(tr3,tr3_,axis=0).astype(int)
        
        times=kwargs.get('times',self._obj.time.values)
        
        z = kwargs.get('z',self._obj[var].values[0,:].flatten())
        
        zmin = self._obj[var].values.min()
        zmax = self._obj[var].values.max()
        
        nodes = pd.DataFrame({'longitude':x,'latitude':y, '{}'.format(var):z})
        elems  = pd.DataFrame(tri3, columns=['a', 'b', 'c'])
                    
        width = kwargs.get('width', 800)
        height = kwargs.get('height', 600)        
        opts.defaults(opts.WMTS(width=width, height=height))

        tiles = gv.WMTS('https://maps.wikimedia.org/osm-intl/{Z}/{X}/{Y}@2x.png')

        points = gv.operation.project_points(gv.Points(nodes, vdims=['{}'.format(var)]))

        trimesh = gv.TriMesh((elems, points))
        
        def time_mesh(time):
            points.data[var] = self._obj[var].sel(time=time).values
            return gv.TriMesh((elems, points))#, crs=ccrs.GOOGLE_MERCATOR)

        meshes = hv.DynamicMap(time_mesh, kdims=['Time']).redim.values(Time=times)
        
        imesh = rasterize(meshes, aggregator='mean').opts(cmap='viridis', colorbar=True, padding=.1, tools=['hover'], clim=(zmin, zmax))
        
        return hv.output(tiles*imesh, holomap='scrubber', fps=1)
        
        
            
     
