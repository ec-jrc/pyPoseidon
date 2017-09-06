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


def anim(x,y,z,title=None,label=None,units=None,**kwargs):
    
    [v1,v2] = kwargs.get('vrange', [None,None])
    
    v = np.linspace(v1, v2, 15, endpoint=True)
    
    def animate(i):
        H.set_array(z[i,:,:])
        return H
    
    parallels=np.arange(-90.,90.,5.)
    meridians=np.arange(0.,360.,5.)
    
    minlon = kwargs.get('lon0', x.min())
    maxlon = kwargs.get('lon1', x.max())
    minlat = kwargs.get('lat0', y.min())
    maxlat = kwargs.get('lat1', y.max())

    latlons = kwargs.get('latlons', True)
    png = kwargs.get('png', False)
    path = kwargs.get('path', './')
    # Create map
    m = Basemap(projection='cyl', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='l')
    
    fig = plt.figure(figsize=(12,8))
    ax=ax = plt.gca()
  
    if latlons:
       m.drawparallels(parallels,labels=[False,True,True,False])
       m.drawmeridians(meridians,labels=[True,False,False,True])
    
    m.drawcoastlines()
    m.drawmapboundary()
    
    # ims is a list of lists, each row is a list of artists to draw in the
    # current frame; here we are animating three artists, the contour and 2 
    # annotatons (title), in each frame
    ims = []
    for i in range(len(z[:,0,0])):
        im = m.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
        add_arts = im.collections
        text = 'time={}'.format(label[i])
        te = ax.text(90, 90, text)
        an = ax.annotate(text, xy=(-0.25, -0.15), xycoords='axes fraction')
        ims.append(add_arts + [te,an])
    
    #cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])    
    cbar = plt.colorbar(ticks=v,orientation='horizontal', extend='both')#,fraction=0.046, pad=0.04)
    if label : cbar.ax.set_xlabel(units)

    if title : plt.title(title)
    
    anim = animation.ArtistAnimation(fig, ims, interval=200, blit=False,repeat=False)
    
#    anim = animation.FuncAnimation(fig,animate,frames=z.shape[0],interval=200, blit=False,repeat=False)
        
    save = kwargs.get('save', False)
    path = kwargs.get('path', './')
    
    if save : 
        anim.save(path, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])
    else:
        return anim
 
 
 
def video(fname, mimetype):
     """Load the video in the file `fname`, with given mimetype, and display as HTML5 video.
     """
     from IPython.display import HTML
     video_encoded = open(fname, "rb").read().encode("base64")
     video_tag = '<video controls alt="test" src="data:video/{0};base64,{1}">'.format(mimetype, video_encoded)
     return HTML(data=video_tag)

 
 
      