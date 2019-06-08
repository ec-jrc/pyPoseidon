"""
Visualization module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import animation
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr

matplotlib.rc('animation',html='html5')
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['animation.embed_limit'] = '200.'
         
def contour(grid_x,grid_y,z,t,**kwargs):
    fig, ax = plt.subplots(figsize=(12,8)) 
    vmin = kwargs.get('vmin', z.min())
    vmax = kwargs.get('vmax', z.max())
    
    nv = kwargs.get('nv', 10)
    
    title = kwargs.get('title', None)
        
    vrange=np.linspace(vmin,vmax,nv,endpoint=True)
    ## CHOOSE YOUR PROJECTION
 #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
    ax = plt.axes(projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.

    ax.set_aspect('equal')
    ims = []
    for i in range(len(t)):
        im = ax.contourf(grid_x, grid_y, z[i,:,:], vrange, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree())
#        im = ax.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
        add_arts = im.collections
        text = 'time={}'.format(t[i])
        #te = ax.text(90, 90, text)
        an = ax.annotate(text, xy=(0.05, 1.05), xycoords='axes fraction')
        ims.append(add_arts + [an])
    if title : ax.set_title(title) 
    #ax.set_global()
    ax.coastlines('50m')
    #ax.set_extent([grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()])


#cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])    
    cbar = fig.colorbar(im,ticks=vrange,orientation='vertical', extend='both')#,fraction=0.046, pad=0.04)
#plt.colorbar()

    v = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat=False)
  
    plt.close()
  
    return v


def update_quiver(num, Q, U, V, step):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """

    Q.set_UVC(U[num,::step,::step],V[num,::step,::step])

    return Q,

def update_qframes(num, Q, U, V):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """

    Q.set_UVC(U[num,:],V[num,:])

    return Q,
    
    
    
def quiver(X,Y,U,V,t,**kwargs):
    
    fig = plt.figure(figsize=(12,8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    crs = ccrs.PlateCarree()
    ax.set_aspect('equal')

    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'],zorder=0)

    sea_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['water'], zorder=0)

    title = kwargs.get('title', None)

    ax.coastlines('50m')
    ax.add_feature(land_50m)
    ax.add_feature(sea_50m)

    scale = kwargs.get('scale', 1.) # change accordingly to fit your needs
    step = kwargs.get('step', 1) # change accordingly to fit your needs

    Q = ax.quiver(X[::step,::step], Y[::step,::step], U[0,::step,::step], V[0,::step,::step], pivot='mid', color='k', angles='xy', scale_units='xy', scale = scale, transform=crs)

    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    ax.set_title(title) 
    #ax.set_global()

    plt.close()
    # you need to set blit=False, or the first set of arrows never gets
    # cleared on subsequent frames
    v = animation.FuncAnimation(fig, update_quiver, fargs=(Q, U, V, step), frames = range(0,np.size(t)),
                                blit=False, repeat=False)#, interval=1)    
    
    return v
    
 
def video(fname, mimetype):
     """Load the video in the file `fname`, with given mimetype, and display as HTML5 video.
     """
     from IPython.display import HTML
     import base64
     video = open(fname, "rb").read()
     video_encoded = base64.b64encode(video)
     video_tag = '<video controls alt="test" src="data:video/{0};base64,{1}">'.format(mimetype, video_encoded.decode('ascii'))
     return HTML(data=video_tag)


@xr.register_dataset_accessor('splot')

class splot(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj    
        
 
    def contour(self,var,**kwargs):        
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        t = kwargs.get('t',self._obj.time.values)
        
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        tri3 = tri3-1 # adjust for python   
         
        if np.abs(x.min()-x.max()) > 359.:
            # Use Matplotlib for triangulation
            triang = matplotlib.tri.Triangulation(x, y)
            tri3 = triang.triangles
    
        it = kwargs.get('it', None)
        z = self._obj[var].values[it,:].flatten()
        
        fig, ax = plt.subplots(figsize=(12,8)) 
        vmin = kwargs.get('vmin', z.min())
        vmax = kwargs.get('vmax', z.max())
    
        nv = kwargs.get('nv', 10)
    
        title = kwargs.get('title', 'contour plot for {}'.format(var))
        
        vrange=np.linspace(vmin,vmax,nv,endpoint=True)
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            z = np.ma.masked_array(z,mask)
            z = z.filled(fill_value=-99999)
        
        plt.gca().set_aspect('equal')
        
        p = plt.tricontour(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree() )
        cbar = fig.colorbar(p,ticks=vrange,orientation='vertical', extend='both')
        if it:
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=(0.05, -0.1), xycoords='axes fraction')
        
        ax.set_title(title,pad=30) 
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        
        ax.coastlines('50m'); ax.gridlines(draw_labels=True);
        
        return p   
    
    def contourf(self,var,**kwargs):
        
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        t = kwargs.get('t',self._obj.time.values)
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        tri3 = tri3-1 # adjust for python    
    
        if np.abs(x.min()-x.max()) > 359.:
            # Use Matplotlib for triangulation
            triang = matplotlib.tri.Triangulation(x, y)
            tri3 = triang.triangles
    
        it = kwargs.get('it', None)
        
        z = self._obj[var].values[it,:].flatten()
        
        fig, ax = plt.subplots(figsize=(12,8)) 
        vmin = kwargs.get('vmin', z.min())
        vmax = kwargs.get('vmax', z.max())
    
        nv = kwargs.get('nv', 10)
    
        title = kwargs.get('title', 'contourf plot for {}'.format(var))
        
        vrange=np.linspace(vmin,vmax,nv,endpoint=True)
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            z = np.ma.masked_array(z,mask)
            z = z.filled(fill_value=-99999)
        
        
        plt.gca().set_aspect('equal')
        
        p = plt.tricontourf(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, transform=ccrs.PlateCarree() )
        cbar = fig.colorbar(p,ticks=vrange,orientation='vertical', extend='both')
        if it :
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=(0.05, -.1), xycoords='axes fraction')
        
        ax.set_title(title,pad=30) 
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        
        ax.coastlines('50m'); ax.gridlines(draw_labels=True);
        
        return p   
            
    
    def quiver(self,var,**kwargs):
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        t = kwargs.get('t',self._obj.time.values)
        
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        tri3 = tri3-1 # adjust for python    
        
        if np.abs(x.min()-x.max()) > 359.:
            # Use Matplotlib for triangulation
            triang = matplotlib.tri.Triangulation(x, y)
            tri3 = triang.triangles
        
        it = kwargs.get('it', None)
        
        u = self._obj[var].values[it,:,0].flatten()
        v = self._obj[var].values[it,:,1].flatten()
        
        scale = kwargs.get('scale', .1)
        
        fig = plt.figure(figsize=(12,8)) 
        title = kwargs.get('title', 'vector plot for {}'.format(var))
        
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            u = np.ma.masked_array(u,mask)
            v = np.ma.masked_array(v,mask)
            v = v.filled(fill_value=-99999)    
            u = u.filled(fill_value=-99999)
        
        
        plt.gca().set_aspect('equal')
        
        p = plt.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=scale)
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')        
        ax.set_title(title, pad=30) 
        
        if it :
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=(0.05, -.1), xycoords='axes fraction')
        
        
        ax.coastlines('50m'); ax.gridlines(draw_labels=True);
        
        
        return p   
        
    def grid(self,**kwargs):
                  
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        tri3 = tri3-1 # adjust for python    

        if np.abs(x.min()-x.max()) > 359.:
            # Use Matplotlib for triangulation
            triang = matplotlib.tri.Triangulation(x, y)
            tri3 = triang.triangles 
         
        fig, ax = plt.subplots(figsize=(12,8))
        
        ax = plt.axes(projection=ccrs.PlateCarree())
        
        plt.gca().set_aspect('equal')
        
        
        ax.set_aspect('equal')    
        g = plt.triplot(x,y,tri3,'go-', lw=.5, markersize=5)#, transform=ccrs.PlateCarree() )
        
        title = kwargs.get('title', 'Grid plot')
        ax.set_title(title, pad=30)
        
        ax.coastlines('50m'); ax.gridlines(draw_labels=True);
        
        
        return g
    
    
    
    def qframes(self,var,**kwargs):
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        
        u = kwargs.get('u',self._obj[var].values[:,:,0])
        v = kwargs.get('v',self._obj[var].values[:,:,1])
        
        t = kwargs.get('t',self._obj.time.values)
        
        scale = kwargs.get('scale', .1)
        
        fig = plt.figure(figsize=(12,8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        crs = ccrs.PlateCarree()
        ax.set_aspect('equal')

        land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'],zorder=0)

        sea_50m = cfeature.NaturalEarthFeature('physical', 'ocean', '50m',
                                                edgecolor='face',
                                                facecolor=cfeature.COLORS['water'], zorder=0)

        title = kwargs.get('title', None)

        ax.coastlines('50m')
        ax.add_feature(land_50m)
        ax.add_feature(sea_50m)

        scale = kwargs.get('scale', 1.) # change accordingly to fit your needs
        step = kwargs.get('step', 1) # change accordingly to fit your needs

        Q = ax.quiver(x, y, u[0,:], v[0,:], pivot='mid', color='k', angles='xy', scale_units='xy', scale = scale, transform=crs)

        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_title(title) 
        #ax.set_global()

        plt.close()
        # you need to set blit=False, or the first set of arrows never gets
        # cleared on subsequent frames
        v = animation.FuncAnimation(fig, update_qframes, fargs=(Q, u, v), blit=False, repeat=False,
                               frames = range(0,np.size(t)))   
        
        return v
                
 
    def frames(self,var,**kwargs):
    
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        t = kwargs.get('t',self._obj.time.values)
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        tri3 = tri3-1 # adjust for python    

        if np.abs(x.min()-x.max()) > 359.:
            # Use Matplotlib for triangulation
            triang = matplotlib.tri.Triangulation(x, y)
            tri3 = triang.triangles
    
        z = self._obj[var].values
        
        fig, ax = plt.subplots(figsize=(12,8)) 
        vmin = kwargs.get('vmin', z.min())
        vmax = kwargs.get('vmax', z.max())
    
        nv = kwargs.get('nv', 10)
    
        title = kwargs.get('title', None)
        
        vrange=np.linspace(vmin,vmax,nv,endpoint=True)
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            z = np.ma.masked_array(z,mask)
            z = z.filled(fill_value=-99999)
        
        
    ## CHOOSE YOUR PROJECTION
 #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
    #ax = plt.axes(projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.

        ax.set_aspect('equal')
        ims = []
        for i in range(len(t)):
            im = ax.tricontourf(x, y, tri3, z[i,:], vrange, vmin=vmin, vmax=vmax)#, transform=ccrs.PlateCarree())
#        im = ax.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
            add_arts = im.collections
            text = 'time={}'.format(t[i])
        #te = ax.text(90, 90, text)
            an = ax.annotate(text, xy=(0.05, -.1), xycoords='axes fraction')
            ims.append(add_arts + [an])
        if title : ax.set_title(title) 
    #ax.set_global()
    #ax.coastlines('50m')
    #ax.set_extent([grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()])


#cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])    
        cbar = fig.colorbar(im,ticks=vrange,orientation='vertical', extend='both')#,fraction=0.046, pad=0.04)
#plt.colorbar()

        v = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat=False)
     
        plt.close()
    
        return v
 
      
