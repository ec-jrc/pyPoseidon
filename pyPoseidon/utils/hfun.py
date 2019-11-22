import pandas as pd
import shapely
import numpy as np
import xarray as xr
from .limgrad import *
from matplotlib import tri


def hfun(dem, path='.', tag='jigsaw', resolution_min=.05, resolution_max=.5, dhdx=.15, imax=100, **kwargs):
    
      
    X, Y = np.meshgrid(dem.Dataset.longitude.values,dem.Dataset.latitude.values)
    V = dem.Dataset.values
    
    #scale
    hmin = resolution_min                       # min. H(X) [deg.]
    hmax = resolution_max                       # max. H(X)

    V[V>0] = 0 #normalize to only negative values

    hfun =  np.sqrt(-V)/.5 # scale with sqrt(H)
    hfun[hfun < hmin] = hmin
    hfun[hfun > hmax] = hmax
    
    hfun = hfun.flatten() # make it 1-d
    
    hfun = np.array([[hf] for hf in list(hfun)]) #convert it to the appropriate format for LIMHFN2 below
    
    #triangulate
    points = np.column_stack([X.flatten(),Y.flatten()])
    # Use Matplotlib for triangulation
    triang = tri.Triangulation(points[:,0], points[:,1])
#    tri3 = triang.triangles
    edges = triang.edges
    #edge lengths
    elen = [shapely.geometry.LineString(points[edge]).length for edge in edges]
    
    [fun,flag] = limgrad2(edges,elen,hfun,dhdx,imax)

    cfun = fun.flatten().reshape(X.shape).T

    ##OUTPUT
    
    dh = xr.Dataset({'z': (['longitude', 'latitude'], cfun)},
                coords={'longitude': ('longitude', dem.Dataset.longitude.values),
                        'latitude': ('latitude', dem.Dataset.latitude.values)})
    
    
    return dh
             
