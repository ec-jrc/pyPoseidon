import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid


def map(x, y, z, title=None,label=None,**kwargs):
    
    minlon = kwargs.get('lon0', x.min())
    maxlon = kwargs.get('lon1', x.max())
    minlat = kwargs.get('lat0', y.min())
    maxlat = kwargs.get('lat1', y.max())

# Create map
    m = Basemap(projection='cyl', llcrnrlat=minlat,urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon,resolution='l')
    fig = plt.figure(figsize=(10,8))
    cs = m.contourf(x,y,z,cmap=plt.cm.jet)
    m.drawcoastlines()
    m.drawmapboundary()
    plt.title(title)
    cbar = plt.colorbar(orientation='horizontal', extend='both')
    cbar.ax.set_xlabel(label)

# Save figure (without 'white' borders)
#plt.savefig('topo.png', bbox_inches='tight')
    plt.show(block=False)
