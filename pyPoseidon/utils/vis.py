"""
Visualization module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
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
import geopandas as gp
import shapely

matplotlib.rc('animation',html='html5')
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['animation.embed_limit'] = '200.'
plt.style.use(['dark_background'])
         
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


@xr.register_dataset_accessor('pplot')
#@xr.register_dataarray_accessor('pplot')

class pplot(object):
    
    def __init__(self, xarray_obj):
        self._obj = xarray_obj    
        
 
    def contour(self,**kwargs):        
        
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
                
        fig, ax = plt.subplots(figsize=(12,8)) 
        vmin = kwargs.get('vmin', z.min())
        vmax = kwargs.get('vmax', z.max())
    
        nv = kwargs.get('nv', 10)
        xy=kwargs.get('xy',(.3,1.05))
        title = kwargs.get('title', 'contour plot for {}'.format(var))
        
        vrange=np.linspace(vmin,vmax,nv,endpoint=True)
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(x.mean(), y.mean()))
    #   ax = plt.axes(projection=ccrs.PlateCarree())
    #   ax.background_patch.set_facecolor('k')
        
        ax = plt.axes()                       
        
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            z = np.ma.masked_array(z,mask)
            z = z.filled(fill_value=-99999)
        
        
        for val in ['x','y','t','it','vmin','vmax','title','nv','tri3', 'mask','xy','z', 'var']:
            try:
                del kwargs[val]
            except:
                pass        
          
        ax.set_aspect('equal')
                
        p = plt.tricontour(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, **kwargs)
        cbar = fig.colorbar(p,ticks=vrange,orientation='vertical', extend='both')
        if it:
    
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=xy, xycoords='axes fraction')
        
        ax.set_title(title,pad=30) 
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
                
        return p, ax   
    
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
    
        nv = kwargs.get('nv', 10)
    
        title = kwargs.get('title', 'contourf plot for {}'.format(var))
        
        vrange=np.linspace(vmin,vmax,nv,endpoint=True)
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))        
    #    [fig,ax] = kwargs.get('figure',[plt.figure(figsize=(12,8)),plt.axes(projection=ccrs.PlateCarree())])
    #    ax.set_extent([x.min(), x.max(), y.min(), y.max()])    
    #     ax.background_patch.set_facecolor('k')
    
        fig = plt.figure(figsize=(12,8))
        ax = plt.axes()                       
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            z = np.ma.masked_array(z,mask)
            z = z.filled(fill_value=-99999)
        
        xy=kwargs.get('xy',(.3,1.05))
        
        for val in ['x','y','t','it','z','vmin','vmax','title','nv','tri3', 'mask','xy','var','figure']:
            try:
                del kwargs[val]
            except:
                pass        
        
        ax.set_aspect('equal')
        
        
        p = ax.tricontourf(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, **kwargs)#, transform=ccrs.PlateCarree() )
        cbar = fig.colorbar(p,ticks=vrange,orientation='vertical', extend='both')
        if it :
                       
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=xy, xycoords='axes fraction')
        
        ax.set_title(title,pad=30) 
        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
        
        return fig, ax   
            
    
    def quiver(self,**kwargs):
                
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        try:
            t = kwargs.get('t',self._obj.time.values)
        except:
            pass
        
        
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
                        
        it = kwargs.get('it', None)
        
        var = kwargs.get('var','dahv')
        u = kwargs.get('u',self._obj[var].values[it,:,0].flatten())
        v = kwargs.get('v',self._obj[var].values[it,:,1].flatten())
        
        scale = kwargs.get('scale', .1)
        color = kwargs.get('color', 'white')
        
        fig = plt.figure(figsize=(12,8)) 
        title = kwargs.get('title', 'vector plot for {}'.format(var))
        xy=kwargs.get('xy',(0.05, -.1))
        
       ## CHOOSE YOUR PROJECTION
    #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
    #   ax = plt.axes(projection=ccrs.PlateCarree())
    #   ax.background_patch.set_facecolor('k')
        
        ax = plt.gca()                       
        
        
        #optional mask for the data
        mask = kwargs.get('mask',None)
        if 'mask' in kwargs:
            u = np.ma.masked_array(u,mask)
            v = np.ma.masked_array(v,mask)
            v = v.filled(fill_value=-99999)    
            u = u.filled(fill_value=-99999)
        
        for val in ['x','y','t','it','u','v','title','tri3', 'xy', 'scale','mask','color','var']:
            try:
                del kwargs[val]
            except:
                pass        
        
        ax.set_aspect('equal')
                
        
        p = plt.quiver(x, y, u, v, angles='xy', scale_units='xy', scale=scale, color=color, **kwargs)
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')        
        ax.set_title(title, pad=30) 
        
        if it :
                     
            text = 'time={}'.format(t[it])
            an = ax.annotate(text, xy=xy, xycoords='axes fraction')
        
        
        return p, ax  
        
    def grid(self, **kwargs):
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
        
        for val in ['x','y','tri3']:
            try:
                del kwargs[val]
            except:
                pass
                
        fig = plt.figure(figsize=(12,8))
        ax = plt.gca()       
        #ax = plt.axes(projection=ccrs.PlateCarree())        
       # ax.background_patch.set_facecolor('k')

        ax.set_aspect('equal')
               
        g = plt.triplot(x,y,tri3,'go-', **kwargs)#, lw=.5, markersize=5)#, transform=ccrs.PlateCarree() )
        
        title = kwargs.get('title', 'Grid plot')
        ax.set_title(title, pad=30)
        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
                
        
        return fig , ax
    
    
    
    def qframes(self,**kwargs):
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        
        
        cr = kwargs.get('coastlines', None)
        c_attrs = kwargs.get('coastlines_attrs', {})
        
        var = kwargs.get('var','dahv')
        u = kwargs.get('u',self._obj[var].values[:,:,0])
        v = kwargs.get('v',self._obj[var].values[:,:,1])
        
        t = kwargs.get('t',self._obj.time.values)
        
        color = kwargs.get('color', 'white')
        
        
#        ax = plt.axes(projection=ccrs.PlateCarree())
      #  ax.set_extent([x.min(), x.max(), y.min(), y.max()])
         
        fig = plt.figure(figsize=(12,8))
        ax = plt.gca()                       
        

        ax.set_aspect('equal')
        
        title = kwargs.get('title', None)

        scale = kwargs.get('scale', 1.) # change accordingly to fit your needs
        step = kwargs.get('step', 1) # change accordingly to fit your needs

        Q = ax.quiver(x, y, u[0,:], v[0,:], pivot='mid', color=color, angles='xy', scale_units='xy', scale = scale)

#        if cr is not None:
#            try:           
#                coastl = gp.GeoDataFrame.from_file(cr)
#            except:
#                coastl = gp.GeoDataFrame(cr)
#            coastl.plot(ax=ax, **c_attrs)
        

        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_title(title) 
        #ax.set_global()

        # you need to set blit=False, or the first set of arrows never gets
        # cleared on subsequent frames
        v = animation.FuncAnimation(fig, update_qframes, fargs=(Q, u, v), blit=False, repeat=False,
                               frames = range(0,np.size(t)))   
        
        plt.close()
        
        return v
                
 
    def frames(self,**kwargs):
    
        cr = kwargs.get('coastlines', None)
        c_attrs = kwargs.get('coastlines_attrs', {})
        
        x = kwargs.get('x',self._obj.SCHISM_hgrid_node_x[:].values)
        y = kwargs.get('y',self._obj.SCHISM_hgrid_node_y[:].values)
        t = kwargs.get('t',self._obj.time.values)
        tri3 = kwargs.get('tri3',self._obj.SCHISM_hgrid_face_nodes.values[:,:3].astype(int))
    
        var = kwargs.get('var','depth')
        z = kwargs.get('z',self._obj[var].values)
        
        # set figure size
        xr = x.max() - x.min() 
        yr = y.max() - y.min()
        ratio = yr/xr
        xf=12
        yf=np.ceil(12*ratio).astype(int)
        
        fig = plt.figure(figsize=(xf,yf)) 
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
      #  ax = plt.axes(projection=ccrs.PlateCarree())
     #   ax.background_patch.set_facecolor('k')        
    # Limit the extent of the map to a small longitude/latitude range.
    #    ax.set_extent([x.min(), x.max(), y.min(), y.max()])
        
        ax = plt.axes()                       
        ax.set_aspect('equal')
        
        ims = []
        for i in range(len(t)):
            im = ax.tricontourf(x, y, tri3, z[i,:], vrange, vmin=vmin, vmax=vmax)#, transform=ccrs.PlateCarree())
#        im = ax.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
            add_arts = im.collections
            text = 'time={}'.format(t[i])
        #te = ax.text(90, 90, text)
            an = ax.annotate(text, xy=(0.05, -.2), xycoords='axes fraction')
            ims.append(add_arts + [an])
            
#            if cr is not None: TO DO
#                try:
#                    coastl = gp.GeoDataFrame.from_file(cr)
#                except:
#                    coastl = gp.GeoDataFrame(cr)
#                coastl.plot(ax=ax, **c_attrs)
                    
        if title : ax.set_title(title) 
    #ax.set_global()
    #ax.coastlines('50m')
        


#cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])    
        cbar = fig.colorbar(im,ticks=vrange,orientation='vertical', extend='both',fraction=0.017, pad=0.04)
#plt.colorbar()

        v = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat=False)
     
        plt.close()
    
        return v
 
