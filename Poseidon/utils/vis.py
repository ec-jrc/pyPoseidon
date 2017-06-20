import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib import animation,rc

rc('animation',html='html5')


def map(x, y, z, title=None,label=None,**kwargs):
    
    parallels=np.arange(-90.,90.,5.)
    meridians=np.arange(0.,360.,5.)
    
    minlon = kwargs.get('lon0', x.min())
    maxlon = kwargs.get('lon1', x.max())
    minlat = kwargs.get('lat0', y.min())
    maxlat = kwargs.get('lat1', y.max())

    ticks = kwargs.get('ticks', True)
    png = kwargs.get('png', False)
    path = kwargs.get('path', './')
# Create map
    m = Basemap(projection='cyl', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='l')
    fig = plt.figure(figsize=(10,8))
    cs = m.contourf(x,y,z,cmap=plt.cm.jet)
    if ticks:
       m.drawparallels(parallels,labels=[False,True,True,False])
       m.drawmeridians(meridians,labels=[True,False,False,True])
    
    m.drawcoastlines()
    m.drawmapboundary()
    plt.title(title)
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel(label)

# Save figure (without 'white' borders)
    if png : plt.savefig(path, bbox_inches='tight')
    plt.show(block=False)


def anim(x,y,z,title=None,label=None,**kwargs):
    
    def animate(i):
        H.set_array(z[i,:,:])
        return H
    
    parallels=np.arange(-90.,90.,5.)
    meridians=np.arange(0.,360.,5.)
    
    minlon = kwargs.get('lon0', x.min())
    maxlon = kwargs.get('lon1', x.max())
    minlat = kwargs.get('lat0', y.min())
    maxlat = kwargs.get('lat1', y.max())

    ticks = kwargs.get('ticks', True)
    png = kwargs.get('png', False)
    path = kwargs.get('path', './')
    # Create map
    m = Basemap(projection='cyl', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='l')
    
    fig = plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111)
    H=m.contourf(x,y,z[0,:,:])
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel(label)
    
    if ticks:
       m.drawparallels(parallels,labels=[False,True,True,False])
       m.drawmeridians(meridians,labels=[True,False,False,True])
    
    m.drawcoastlines()
    m.drawmapboundary()
    
    anim = animation.FuncAnimation(fig,animate,frames=z.shape[0],interval=200, blit=False,repeat=False)
        
    save = kwargs.get('save', False)
    path = kwargs.get('path', './')
    
    if save : anim.save(path, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])
    
    return anim
    