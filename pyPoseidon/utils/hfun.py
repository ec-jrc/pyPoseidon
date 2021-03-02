import pandas as pd
import geopandas as gp
import shapely
import numpy as np
import xarray as xr
from .limgrad import *
import matplotlib
from pyPoseidon.utils.stereo import to_lat_lon, to_stereo
import pyPoseidon
import math

#https://stackoverflow.com/questions/44934631/making-grid-triangular-mesh-quickly-with-numpy
def MakeFacesVectorized1(Nr,Nc):

    out = np.empty((Nr-1,Nc-1,2,3),dtype=int)

    r = np.arange(Nr*Nc).reshape(Nr,Nc)

    out[:,:, 0,0] = r[:-1,:-1]
    out[:,:, 1,0] = r[:-1,1:]
    out[:,:, 0,1] = r[:-1,1:]

    out[:,:, 1,1] = r[1:,1:]
    out[:,:, :,2] = r[1:,:-1,None]

    out.shape =(-1,3)
    return out


def hfun(data, path='.', tag='jigsaw', resolution_min=.05, resolution_max=.5, dhdx=.15, imax=100, **kwargs):


    X, Y = np.meshgrid(data.longitude.values,data.latitude.values)
    V = data.values

    #scale
    hmin = resolution_min                       # min. H(X) [deg.]
    hmax = resolution_max                       # max. H(X)

    V[V>0] = 0 #normalize to only negative values

    hfun =  np.sqrt(-V)/.5 # scale with sqrt(H)
    hfun[hfun < hmin] = hmin
    hfun[hfun > hmax] = hmax

#    subspace = kwargs.get('subspace', None)

#    if subspace is not None:
#        mask = ...
#        hfun[mask] =

    hfun = hfun.flatten() # make it 1-d

    hfun = hfun.reshape(hfun.shape[0],-1) #convert it to the appropriate format for LIMHFN2 below

    #triangulate
    points = np.column_stack([X.flatten(),Y.flatten()])
    tria = MakeFacesVectorized1(V.shape[0],V.shape[1])
    # Use Matplotlib for triangulation
    triang = matplotlib.tri.Triangulation(points[:,0], points[:,1], tria)
#    tri3 = triang.triangles
    edges = triang.edges
    #edge lengths
    ptdiff = lambda p: (p[0][0]-p[1][0], p[0][1]-p[1][1])
    diffs = map(ptdiff , points[edges])
    elen = [math.hypot(d1,d2) for d1,d2 in diffs]

    [fun,flag] = limgrad2(edges,elen,hfun,dhdx,imax)

    cfun = fun.flatten().reshape(X.shape).T

    ##OUTPUT

    dh = xr.Dataset({'z': (['longitude', 'latitude'], cfun)},
                coords={'longitude': ('longitude', data.longitude.values),
                        'latitude': ('latitude', data.latitude.values)})


    return dh


def to_hfun_mesh(dh,fhfun):

    dps = dh[['u','v','z','h']].to_dataframe()
    tria3 = dh.tria.to_pandas()

    with open(fhfun,'w') as f:
        f.write('#{}; created by pyPoseidon\n'.format(pyPoseidon.__version__))
        f.write('MSHID=3;EUCLIDEAN-MESH\n')
        f.write('NDIMS=2\n')
        f.write('POINT={}\n'.format(dps.shape[0]))

    with open(fhfun, 'a') as f:
        dps[['u','v','z']].to_csv(f, index=False, header=0, sep=';')

    with open(fhfun, 'a') as f:
        f.write('VALUE={};1\n'.format(dps.shape[0]))
        dps[['h']].to_csv(f, index=False, header=0)


    with open(fhfun, 'a') as f:
        f.write('TRIA3={}\n'.format(tria3.shape[0]))
        tria3.to_csv(f, index=False, header=0, sep=';')



def to_hfun_grid(dh,fhfun):
        # write hfun file

        # write header
        with open(fhfun,'w') as f:
            f.write('#{}; created by pyPoseidon\n'.format(pyPoseidon.__version__))
            f.write('MSHID=3;EUCLIDEAN-GRID\n')
            f.write('NDIMS=2\n')
            f.write('COORD=1;{}\n'.format(dh.longitude.size))

        with open(fhfun, 'a') as f:
            np.savetxt(f, dh.longitude.values)

        with open(fhfun, 'a') as f:
            f.write('COORD=2;{}\n'.format(dh.latitude.size))

        with open(fhfun, 'a') as f:
            np.savetxt(f, dh.latitude.values)

        with open(fhfun, 'a') as f:
            f.write('VALUE={};1\n'.format(dh.z.size))

        with open(fhfun, 'a') as f:
            np.savetxt(f, dh.z.values.flatten())



def hfun_(coastlines,res=.1, R=1.):

    amask = (coastlines.bounds.miny < coastlines.bounds.miny.min() + .1)
    anta = coastlines[amask]
    anta = anta.reset_index(drop=True)

    ### convert to stereo
    try:
        ant = pd.DataFrame(anta.boundary.values[0].coords[:], columns=['lon','lat'])
    except:
        ant = pd.DataFrame(anta.boundary.explode().values[0].coords[:], columns=['lon','lat']) # convert boundary values to pandas

    d1 = ant.where(ant.lon==ant.lon.max()).dropna().index[1:] # get artificial boundaries as -180/180
    d2 = ant.where(ant.lon==ant.lon.min()).dropna().index[1:]
    ant = ant.drop(d1).drop(d2) # drop the points
    d3 = ant.where(ant.lat==ant.lat.min()).dropna().index # drop lat=-90 line
    ant = ant.drop(d3)
    ub, vb = to_stereo(ant.lon.values,ant.lat.values, R)
    ant.lon=ub
    ant.lat=vb

    an = gp.GeoDataFrame({'geometry' : [shapely.geometry.LineString(ant.values)], 'length':shapely.geometry.LineString(ant.values).length }) # put together a LineString

    # create simple grid
    d1 = np.linspace(an.bounds.minx,an.bounds.maxx,100)
    d2 = np.linspace(an.bounds.miny,an.bounds.maxy,100)

    ui, vi = np.meshgrid(d1,d2)
    # Use Matplotlib for triangulation
    triang = matplotlib.tri.Triangulation(ui.flatten(), vi.flatten())
    tri = triang.triangles


    #stereo->2D scale
    ci=4*R**2/(ui**2+vi**2+4*R**2)
    ci

    # create weight field
    points = np.column_stack([ui.flatten(),vi.flatten()])
    dps = pd.DataFrame(points,columns=['u','v'])
    dps['z'] = 0
    dps['h'] = res / ci.flatten()

    tria3 = pd.DataFrame(tri,columns=['a','b','c'])
    tria3['w'] = 0

    p1 = dps.to_xarray()
    p1 = p1.rename({'index':'nodes'})
    p1 = p1.assign({'tria':(['elem','n'], tria3.values)})

    return p1


