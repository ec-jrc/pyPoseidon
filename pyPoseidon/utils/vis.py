import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib import animation,rc
from cartopy import crs
import cartopy.crs as ccrs
import cartopy.feature as cfeature


rc('animation',html='html5')
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams['animation.embed_limit'] = '200.'
        
FFWriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264','-pix_fmt','yuv420p'])
 
def contour(grid_x,grid_y,z,t,**kwargs):
    fig, ax = plt.subplots(figsize=(12,8)) 
    [v1,v2] = kwargs.get('vrange', [None,None])
    
    title = kwargs.get('title', None)
    
    if not v1:
        v1=z.min()
        v2=z.max()
    
    v=np.linspace(z.min(),z.max(),10,endpoint=True)
    ## CHOOSE YOUR PROJECTION
 #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
    ax = plt.axes(projection=ccrs.PlateCarree())
    # Limit the extent of the map to a small longitude/latitude range.

    ax.set_aspect('equal')
    ims = []
    for i in range(len(t)):
        im = ax.contourf(grid_x, grid_y, z[i,:,:], v, vmin=v1, vmax=v2, transform=ccrs.PlateCarree())
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
    cbar = fig.colorbar(im,ticks=v,orientation='vertical', extend='both')#,fraction=0.046, pad=0.04)
#plt.colorbar()

    v = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat=False)
     
    plt.close()
    
    if 'savepath' in kwargs.keys(): 
        savepath = kwargs.get('path', './')
        v.save(path, writer = FFWriter)      
    else:
        return v


def update_quiver(num, Q, U, V, step):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """
    Q.set_UVC(U[num,::step,::step],V[num,::step,::step])

    return Q,
    
    
def quiver(X,Y,U,V,t,**kwargs):
    
    fig = plt.figure(figsize=(16,12))
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
    return animation.FuncAnimation(fig, update_quiver, fargs=(Q, U, V, step),
                               interval=np.size(t), blit=False, repeat=False)    
    
 
def video(fname, mimetype):
     """Load the video in the file `fname`, with given mimetype, and display as HTML5 video.
     """
     from IPython.display import HTML
     video_encoded = open(fname, "rb").read().encode("base64")
     video_tag = '<video controls alt="test" src="data:video/{0};base64,{1}">'.format(mimetype, video_encoded)
     return HTML(data=video_tag)

 
 
      
