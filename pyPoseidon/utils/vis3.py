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
                                    
        mlab.show()
        return
                
