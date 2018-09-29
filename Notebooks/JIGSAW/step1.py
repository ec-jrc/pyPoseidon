#!/Users/brey/miniconda3/envs/pyPoseidon/bin/python
# coding: utf-8
"""

Jigsaw Global Grid

This is the prototype for creating a bathymetry weighted global grid.

Step 1. Create the Global grid (full globe coverage)

"""
# #  

# This is the prototype for creating a bathymery weighted global grid.

import pandas as pd
import numpy as np
from pyPoseidon.dem import erdap, gebco
import os
import subprocess


PATH = '/Users/brey/GLOBAL/JIGSAW/'

# ## geo file

fgeo=PATH + 'geo/Global-GEOM.msh'

# write header
with open(fgeo,'w') as f:
    f.write('#{}; created by pyPoseidon\n'.format(fgeo.split('/')[-1]))
    f.write('MSHID=3;ELLIPSOID-MESH\n')
    f.write('RADII={};{};{}\n'.format(6371.000000,6371.000000,6371.000000))


# hfun file


fbath = '/Users/brey/GitHub/jigsaw-geo-matlab/jigsaw/geo/topo.msh'
coord, nc = pd.read_csv(fbath,header=2,skiprows=1,nrows=0,delimiter='=')
nc0, nc1 = nc.split(';')
nlon = int(nc1)
#read lon for all boundary nodes
xb = pd.read_csv(fbath,skiprows=4,header=None,nrows=nlon,names=['lon'])
coord, nc = pd.read_csv(fbath,header=4,skiprows=nlon,nrows=0,delimiter='=')
nc0, nc2 = nc.split(';')
nlat = int(nc2)

#read lat,depth for boundary nodes
yb = pd.read_csv(fbath,skiprows=5 + nlon,header=None,nrows=nlat,names=['lat'])
value, nc = pd.read_csv(fbath,header=4,skiprows=nlon+nlat+1,nrows=0,delimiter='=')

val, val0 = nc.split(';')
val = int(val)
#read 
bath = pd.read_csv(fbath,skiprows=6 + nlon + nlat,header=None,nrows=val,names=['value'])

vb = bath.values.reshape(nlon,nlat)
zlev=vb.T # to match the matlab configuration in which the file fbath was saved on
X, Y = np.meshgrid(xb,yb)

eps = np.finfo(float).eps

# SET RESOLUTION
radE = 6371.0E+00 
      
hfn0 = +50. ;                      # global spacing
hfn2 = +0.;                        # adapt. spacing

 
dhdx = +.10;                        # max. gradients    

hfun = hfn0*np.ones([nlat,nlon]) # note the shape!


b=zlev.copy() # make a copy so that we don't change zlev

b[b>0] = eps #normalize to only negative values

htop = np.sqrt(-b)/1.5 # scale with sqrt(H)
htop[htop<hfn2] = hfn2 # normalize



## optional focused area
#hfn3 = +50.;                        # arctic spacing
#htop[htop>hfn3] = hfn3


htop[zlev>0.] = hfn0 


#hfun[Y>+50.] = htop[Y>+50.] # with focused area
hfun = htop # global

# ## SPHGRID

ALON = X.T * np.pi / 180. 
ALAT = Y.T * np.pi / 180. 
    
#----------------------------------------- spheroidal coord.
xpos = radE * np.cos(ALON) * np.cos(ALAT) 
ypos = radE * np.sin(ALON) * np.cos(ALAT) 
zpos = radE * np.sin(ALAT)  


ppos = pd.DataFrame({'a':xpos.flatten(),'b':ypos.flatten(), 'c':zpos.flatten()})

grid = np.arange(nlat*nlon).reshape(nlon,nlat)


grid = grid.T

quad = np.zeros([(nlat-1)*(nlon-0),4])

nex = 0
for jpos in range(1,nlat):
       
            quad[nex,0] = grid[jpos-1,nlon-1]
            quad[nex,1] = grid[jpos-1,0]
            quad[nex,2] = grid[jpos-0,0]
            quad[nex,3] = grid[jpos-0,nlon-1]
            
            nex = nex+1
            

        
for ipos in range(1,nlon):
    for jpos in range(1,nlat):
       
            quad[nex,0] = grid[jpos-1,ipos-1]
            quad[nex,1] = grid[jpos-1,ipos-0]
            quad[nex,2] = grid[jpos-0,ipos-0]
            quad[nex,3] = grid[jpos-0,ipos-1]
            
            nex = nex+1 
     
quad = quad.astype(int)


# ### LIMHFUN

edge = np.vstack([quad[:,[0,1]],
        quad[:,[1,2]],
        quad[:,[2,3]],
        quad[:,[3,0]]
           ])


edge.sort(axis=1)

edge = np.unique(edge,axis=0)

evec = ppos.values[edge[:,1],:] - ppos.values[edge[:,0],:] 

elen = np.sqrt(np.sum(evec**2,axis=1))


# ### LIMGRAD


def limgrad(edge,elen,ffun,dfdx,imax):
    
    rfun = ffun.T.flatten()
    
    rfun = np.array([[hf] for hf in list(rfun)])
    
    eps = np.finfo(float).eps
    
    nnod = rfun.size
        
    #-- IVEC(NPTR(II,1):NPTR(II,2)) are edges adj. to II-TH node
    nvec = np.hstack([edge[:,0], edge[:,1]])

    ivec = np.hstack([np.arange(edge.shape[0]),np.arange(edge.shape[0])])


    nvec_ = np.sort(nvec, kind='mergesort')
    pidx = np.argsort(nvec, kind='mergesort') # to match with matlab/octave -> https://stackoverflow.com/questions/39484073/matlab-sort-vs-numpy-argsort-how-to-match-results
    ivec = ivec[pidx]
    nvec = nvec_
    
    
    mark = np.full(rfun.size, False, dtype=bool)
    mark[edge[:,0]]=True
    mark[edge[:,1]]=True

    dif = [nvec[i+1]-nvec[i] for i in range(nvec.shape[0]-1)]
    idxx = np.argwhere(np.array(dif) > 0).flatten()
    
    nptr=np.zeros((mark.size,2))
    nptr[mark,0] = np.append(np.array([0]),idxx+1)
    nptr[mark,1] = np.append(idxx, nnod-1)
    
    nptr = nptr.astype(int)

#----------------------------- ASET=ITER if node is "active"
    aset = np.zeros(nnod)
    
    ftol = min(rfun.flatten()) * np.sqrt(eps)
    

#----------------------------- exhaustive 'til all satisfied 
    
    for i in range(1,imax):
    
    #------------------------- find "active" nodes this pass
        aidx = np.argwhere(aset == i - 1 ) 
        aidx = aidx.flatten()
        
        if not aidx.any(): break
      
    #------------------------- reorder => better convergence

        aval = np.sort(rfun.reshape(ffun.shape).flatten()[aidx],kind='mergesort')
        idxx = np.argsort(rfun.reshape(ffun.shape).flatten()[aidx], kind='mergesort')
        
        aidx = aidx[idxx]
       
    #%------------------------- visit adj. edges and set DFDX
        for ipos in range(len(aidx)):
            npos = aidx[ipos]
        
            for jpos in range(nptr[npos,0], nptr[npos,1]+1):
                
                epos = ivec[jpos]
                
                nod1 = edge[epos,0]
                nod2 = edge[epos,1]
                
               # print ipos, jpos, epos, nod1, nod2

            #----------------- calc. limits about min.-value
                if rfun[nod1] > rfun[nod2]:
                
                    
                    fun1 = rfun[nod2] + elen[epos] * dfdx 
                
                    if rfun[nod1] > fun1+ftol :
                        rfun[nod1] = fun1
                        aset[nod1] = i
                else:
                
                    fun2 = rfun[nod1] + elen[epos] * dfdx 
                    
                    if   rfun[nod2] > fun2+ftol :
                        rfun[nod2] = fun2
                        aset[nod2] = i
     
    flag = i < imax
    
    return rfun,flag


[fun,flag] = limgrad(edge,elen,hfun,dhdx,4096)

tfun = fun.reshape(hfun.T.shape) # we use the transpose to get back to python

#difference from the original hfun
dff = tfun - hfun.T


#MSH FILE

ffun= PATH + 'out/GLOBAL-HFUN.msh'


gf = pd.DataFrame({'hmatx' : xpos.flatten()*np.pi/180, 
'hmaty' : ypos.flatten()*np.pi/180, 
'hmatv' : tfun.flatten()}) 
    

xb['x'] = xb.lon*np.pi/180

yb['y'] = yb.lat*np.pi/180


# write header
with open(ffun,'w') as f:
    f.write('#{}; created by pyPoseidon\n'.format(ffun.split('/')[-1]))
    f.write('MSHID=3;ELLIPSOID-GRID\n')
    f.write('NDIMS=2\n')
    f.write('COORD=1;{}\n'.format(xb.x.size))


with open(ffun, 'a') as f:
    xb.x.to_csv(f, index=False, header=0)


with open(ffun, 'a') as f:
    f.write('COORD=2;{}\n'.format(yb.y.size))


with open(ffun, 'a') as f:
    yb.y.to_csv(f, index=False, header=0)


with open(ffun, 'a') as f:
    f.write('VALUE={};1\n'.format(xb.x.size * yb.y.size)) 


with open(ffun, 'a') as f:
    for i in range(fun.size):
        f.write('{}\n'.format(tfun.flatten()[i]))


# JIG FILE

jig=PATH + 'Global.jig'

fmsh = PATH + 'out/Global.msh'

with open(jig,'w') as f:
    f.write('#{}; created by pyPoseidon\n'.format(jig.split('/')[-1]))
    f.write('GEOM_FILE ={}\n'.format(fgeo))
    f.write('MESH_FILE ={}\n'.format(fmsh))
    f.write('HFUN_FILE ={}\n'.format(ffun))
    f.write('HFUN_SCAL = ABSOLUTE\n')
    f.write('HFUN_HMAX = 150\n')
    f.write('HFUN_HMIN = 0.0\n')
    f.write('MESH_DIMS = 2\n')
    f.write('OPTM_QLIM=0.9375\n')
    f.write('VERBOSITY = 1')


# Execute jigsaw
ex=subprocess.Popen(args=['/Users/brey/GitHub/jigsaw-geo-matlab/jigsaw/bin/MAC-64/jigsaw64r', PATH + 'Global.jig'], cwd='.', shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)

for line in iter(ex.stderr.readline,b''): 
            sys.stdout.write(line)
            sys.stdout.flush()  
ex.stderr.close()            

for line in iter(ex.stdout.readline,b''): 
            sys.stdout.write(line)
            sys.stdout.flush()  
ex.stdout.close()         

