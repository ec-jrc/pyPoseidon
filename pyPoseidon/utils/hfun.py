import pandas as pd
import shapely
import numpy as np

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

    ##OUTPUT
    
    fhfun = path + tag+'-hfun.msh'
    # write header
    with open(fhfun,'w') as f:
        f.write('#{}; created by pyPoseidon\n'.format(tag+'-hfun.msh'))
        f.write('MSHID=3;EUCLIDEAN-GRID\n')
        f.write('NDIMS=2\n')
        f.write('COORD=1;{}\n'.format(dem.Dataset.longitude.size))
        
    with open(fhfun, 'a') as f:
        dem.Dataset.longitude.to_dataframe().to_csv(f, index=False, header=0)
        
    with open(fhfun, 'a') as f:
        f.write('COORD=2;{}\n'.format(dem.Dataset.latitude.size))
        
    with open(fhfun, 'a') as f:
        dem.Dataset.latitude.to_dataframe().to_csv(f, index=False, header=0)
    
    with open(fhfun, 'a') as f:
        f.write('VALUE={};1\n'.format(dem.Dataset.longitude.size * dem.Dataset.latitude.size))
    
    #converted TRANSPOSED, IS THERE A WAY TO FIX IT?
    cfun = fun.flatten().reshape(X.shape).T
    with open(fhfun, 'a') as f:
        for i in range(fun.size):
            f.write('{}\n'.format(cfun.flatten()[i]))
            
