"""
Visualization module in 3D

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

from mayavi import mlab
from mayavi.sources.builtin_surface import BuiltinSurface
import xarray as xr
import numpy as np
import pandas as pd


@xr.register_dataset_accessor('pplot3')
#@xr.register_dataarray_accessor('pplot')

class pplot3(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj    
        

    def globe(self,R):
    # We use a sphere Glyph, throught the points3d mlab function, rather than
    # building the mesh ourselves, because it gives a better transparent
    # rendering.
        sphere = mlab.points3d(0, 0, 0, scale_mode='none',
                                    scale_factor=2*R,
    #                                color=(0.67, 0.77, 0.93),
                                    color=(0, 0., 0.),
                                    resolution=50,
                                    opacity=1.0,
                                    name='Earth')

    # These parameters, as well as the color, where tweaked through the GUI,
    # with the record mode to produce lines of code usable in a script.
        sphere.actor.property.specular = 0.45
        sphere.actor.property.specular_power = 5
    # Backface culling is necessary for more a beautiful transparent
    # rendering.
        sphere.actor.property.backface_culling = True

        return sphere
    
      
    def contourf(self,**kwargs):
                        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        try:
            t = kwargs.get('t',self._obj.time.values)
        except:
            pass
        
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
              
        it = kwargs.get('it', None)
      
        var = kwargs.get('var','depth')
        z = kwargs.get('z',self._obj[var].values[it,:].flatten())
                
        vmin = kwargs.get('vmin', z.min())
        vmax = kwargs.get('vmax', z.max())
        
        R = kwargs.get('R',1.)
        
        px=np.cos(y/180*np.pi)*np.cos(x/180*np.pi)*R
        py=np.cos(y/180*np.pi)*np.sin(x/180*np.pi)*R
        pz=np.sin(y/180*np.pi)*R
        
        rep=kwargs.get('representation','surface')
        
        cmap = kwargs.get('cmap','gist_earth')
        
        mlab.figure(1, size=(3840, 2160), bgcolor=(0, 0, 0), fgcolor=(1.,1.,1.))
        mlab.clf()
        self.globe(R - .02)
        # 3D triangular mesh surface (like trisurf)
        grd = mlab.triangular_mesh(px,py,pz,tri3, representation=rep, opacity=1.0, scalars=z,  colormap=cmap,vmin=vmin,vmax=vmax)
                            
        grd.actor.mapper.scalar_visibility = True
        mlab.view(azimuth=0, distance=4)
        
        title = kwargs.get('title', '{}'.format(var))
        
        mlab.colorbar(grd, title='Bathymetry', orientation='vertical')
        
        coast = kwargs.get('coastlines',None)
        
        if coast is not None :
            src, lines = self.c3d(coast,R=R)
            mlab.pipeline.surface(src, color=(1,0,0), line_width=10, opacity=0.8)
        
        
        mlab.show()
        return
    
    
    def grid(self,**kwargs):
                        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        try:
            t = kwargs.get('t',self._obj.time.values)
        except:
            pass
        
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
              
        R = kwargs.get('R',1.)
        
        px=np.cos(y/180*np.pi)*np.cos(x/180*np.pi)*R
        py=np.cos(y/180*np.pi)*np.sin(x/180*np.pi)*R
        pz=np.sin(y/180*np.pi)*R
        
        mlab.figure(1, size=(3840, 2160), bgcolor=(0, 0, 0), fgcolor=(1.,1.,1.))
        mlab.clf()
        self.globe(R - .02)
        # 3D triangular mesh surface (like trisurf)
        grd = mlab.triangular_mesh(px,py,pz,tri3, representation='wireframe', opacity=1.0)
        
        coast = kwargs.get('coastlines',None)
        
        if coast is not None :
            src, lines = self.c3d(coast,R=R)
            mlab.pipeline.surface(src, color=(1,0,0), line_width=10, opacity=0.8)
                                    
        mlab.show()
        return
                

    def c3d(self,coastlines,R=1):
        
        bo = coastlines.geometry.values
        
        dic={}
        for l in range(len(bo)):
        #    print(l)
            lon=[]
            lat=[]
            try:
                for x,y in bo[l].boundary.coords[:]: 
                    lon.append(x)
                    lat.append(y)
            except:
                for x,y in bo[l].boundary[0].coords[:]: 
                    lon.append(x)
                    lat.append(y)

        
            dic.update({'line{}'.format(l):{'lon':lon,'lat':lat}})

        dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()}
        dff = pd.concat(dict_of_df, axis=0)
        dff['z']=0
        dff.head()
        
        # add 3D coordinates
        dff['x']=np.cos(dff.lat/180*np.pi)*np.cos(dff.lon/180*np.pi)*R
        dff['y']=np.cos(dff.lat/180*np.pi)*np.sin(dff.lon/180*np.pi)*R
        dff['z']=np.sin(dff.lat/180*np.pi)*R
        
        # We create a list of positions and connections, each describing a line.
        # We will collapse them in one array before plotting.
        x = list()
        y = list()
        z = list()
        s = list()
        connections = list()

        # The index of the current point in the total amount of points
        index = 0

        # Create each line one after the other in a loop
        for key, sdf in dff.groupby(level=0):
            x.append(sdf.x.values)
            y.append(sdf.y.values)
            z.append(sdf.z.values)
            N = sdf.shape[0]
            #s.append(np.linspace(-2 * np.pi, 2 * np.pi, N))
            # This is the tricky part: in a line, each point is connected
            # to the one following it. We have to express this with the indices
            # of the final set of points once all lines have been combined
            # together, this is why we need to keep track of the total number of
            # points already created (index)
            connections.append(np.vstack(
                               [np.arange(index,   index + N - 1.5),
                                np.arange(index + 1, index + N - .5)]
                                    ).T)
            index += N

        # Now collapse all positions, scalars and connections in big arrays
        x = np.hstack(x)
        y = np.hstack(y)
        z = np.hstack(z)
        #s = np.hstack(s)
        connections = np.vstack(connections)

        # Create the points
        src = mlab.pipeline.scalar_scatter(x, y, z)#, s)

        # Connect them
        src.mlab_source.dataset.lines = connections
        src.update()

        # The stripper filter cleans up connected lines
        lines = mlab.pipeline.stripper(src)
    
        return src, lines