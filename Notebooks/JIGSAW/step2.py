#!/Users/brey/miniconda3/envs/pyPoseidon/bin/python
# coding: utf-8
"""

Jigsaw Global Grid

This is the prototype for converting .msh grid to .gr3.

Step 2. Cut out land and adapt for SCHISM

"""

import pandas as pd
import numpy as np
import pyresample
from shutil import copyfile
from shapely import geometry, ops
import geopandas as gp
import matplotlib.path as mpltPath



PATH = '/Users/brey/GLOBAL/JIGSAW/'

# ## read .msh
msh=PATH + 'out/Global.msh'

#extract number of elements, number of nodes
p,nn=pd.read_csv(msh,header=2,skiprows=1,nrows=0,delimiter='=')
nn=int(nn)

#read lon,lat,depth for all nodes
xy=pd.read_csv(msh,skiprows=4,header=None,delimiter=';',engine='python',nrows=nn,names=['x','y','z','q'])
xy.index=xy.index+1 # adapt to index start = 1


typ,ne=pd.read_csv(msh,header=4,skiprows=nn,nrows=0,delimiter='=')
ne = int(ne)

tr=pd.read_csv(msh,skiprows=nn+5,header=None,delimiter=';',engine='python',nrows=ne,names=['a','b','c','q'])
tr.index=tr.index+1

tr['q']=3

tri3 = tr.loc[:,['a','b','c']].values

# -> to lat / lon

rad = xy['x']**2+xy['y']**2+xy['z']**2

eps=np.finfo(np.float).eps

rad = np.sqrt(rad)

lat = np.arcsin(xy['z']/rad)
lon = np.arctan2(xy['y'],xy['x'])

lat = lat * 180 / np.pi
lon = lon * 180 / np.pi

xp=np.cos(lat/180*np.pi)*np.cos(lon/180*np.pi)
yp=np.cos(lat/180*np.pi)*np.sin(lon/180*np.pi)
zp=np.sin(lat/180*np.pi)



#xa = xy['x'].values/rad
#ya = xy['y'].values/rad
#za = xy['z'].values/rad

#pd.DataFrame({'xa':xa,'xp':xp})

xy['lon']=lon
xy['lat']=lat


# FIX Grid

# STEP 1

#make copies
grid = xy.copy()
grid = grid.loc[:,:].astype(float)

c3 = tr.copy()
c3.loc[:,['a']]=c3.loc[:,['a']]+1 # we add one to get the range starting from 1 as the index
c3.loc[:,['b']]=c3.loc[:,['b']]+1
c3.loc[:,['c']]=c3.loc[:,['c']]+1

# Compute orientation

c3['x1'] = grid.loc[c3['a'].values,'lon'].values

c3['y1'] = grid.loc[c3['a'].values,'lat'].values

c3['x2'] = grid.loc[c3['b'].values,'lon'].values

c3['y2'] = grid.loc[c3['b'].values,'lat'].values

c3['x3'] = grid.loc[c3['c'].values,'lon'].values

c3['y3'] = grid.loc[c3['c'].values,'lat'].values

c3['val'] = (c3['y2'] - c3['y1']) * (c3['x3'] - c3['x2']) - (c3['x2'] - c3['x1']) * (c3['y3'] - c3['y2'])

# Fix the orientation 

cc = c3.copy() # make a copy

#swap values 2,3 to change orientation
cc.loc[cc.val > 0, 'b'] = c3.loc[c3.val > 0, 'c'].values.astype(int)

cc.loc[cc.val > 0, 'c'] = c3.loc[c3.val > 0, 'b'].values.astype(int)

cc.loc[:,'b'] = cc.loc[:,'b'].astype(int) # make int


# compute again orientation to verify

cc['x1'] = grid.loc[cc['a'].values,'lon'].values 

cc['y1'] = grid.loc[cc['a'].values,'lat'].values 

cc['x2'] = grid.loc[cc['b'].values,'lon'].values

cc['y2'] = grid.loc[cc['b'].values,'lat'].values

cc['x3'] = grid.loc[cc['c'].values,'lon'].values 

cc['y3'] = grid.loc[cc['c'].values,'lat'].values

cc['val'] = (cc['y2'] - cc['y1']) * (cc['x3'] - cc['x2']) - (cc['x2'] - cc['x1']) * (cc['y3'] - cc['y2'])

print cc.val.min(), cc.val.max() # all negative -> counter-clockwise


# STEP 2

# Because of the international line the orientation is not computed correcty for the elements that cross it. So we need to adjust.

# store the sign of lon
cc['s1']=np.sign(cc.x1)
cc['s2']=np.sign(cc.x2)
cc['s3']=np.sign(cc.x3)

p = cc.loc[(abs(cc.x1 - cc.x2) > 180) | (abs(cc.x2 - cc.x3) > 180) | (abs(cc.x1 - cc.x3) > 180) ]# store values

c2 = cc.copy() #make a copy

# for use below
xn = ['a','b','c'] 
dic={'a':'x1','b':'x2','c':'x3'}

#q=[]
#for slist in p.loc[:,['s1','s2','s3']].values:
#    sw = [xn[k] for k in np.argwhere(slist == slist.sum()).flatten()]
#    q.append(sw)

# Use the sign of lon to find the one node that crosses # there might be a better way
q1=[]
q2=[]
for slist in p.loc[:,['s1','s2','s3']].values:
    sw1 = [xn[k] for k in np.argwhere(slist > 0).flatten()]
    q1.append(sw1)
    sw2 = [xn[k] for k in np.argwhere(slist < 0).flatten()]
    q2.append(sw2)

#reshuffle
l=0
for i in p.index:
    lst = max(q1[l],q2[l],key=len)
    v = [item for item in xn if item not in max(q1[l],q2[l],key=len)][0] 
#    print (v)
#    print (p.loc[i,['s1','s2','s3']].sum())
    c2.loc[i, dic[v] ] = c2.loc[i, dic[v] ] + p.loc[i,['s1','s2','s3']].sum() * 360.
    l+=1

c2['val'] = (c2['y2'] - c2['y1']) * (c2['x3'] - c2['x2']) - (c2['x2'] - c2['x1']) * (c2['y3'] - c2['y2'])

c2['s1']=np.sign(c2.x1)

c2['s2']=np.sign(c2.x2)

c2['s3']=np.sign(c2.x3)

c22 = c2.copy()

#swap values 2,3 to change orientation
c22.loc[c22.val > 0, 'b'] = cc.loc[c2.val > 0, 'c'].values.astype(int)

c22.loc[c22.val > 0, 'c'] = cc.loc[c2.val > 0, 'b'].values.astype(int)

c22.loc[c22.val > 0, 'x2'] = cc.loc[c2.val > 0, 'x3'].values

c22.loc[c22.val > 0, 'x3'] = cc.loc[c2.val > 0, 'x2'].values

c22.loc[c22.val > 0, 'y2'] = cc.loc[c2.val > 0, 'y3'].values

c22.loc[c22.val > 0, 'y3'] = cc.loc[c2.val > 0, 'y2'].values



c22['s1']=np.sign(c22.x1)

c22['s2']=np.sign(c22.x2)

c22['s3']=np.sign(c22.x3)


c22['val'] = (c22['y2'] - c22['y1']) * (c22['x3'] - c22['x2']) - (c22['x2'] - c22['x1']) * (c22['y3'] - c22['y2'])


pf = c22.loc[c22.val > 0]


# ## one more issue (sometimes)

r = pf.loc[pf.val == pf.loc[pf.val < 1000.,'val'].min()] # there is something about node 1 CHECK

#swap again to fix it for good ?
c22.loc[r.index, 'b'] = c2.loc[r.index, 'b'].values.astype(int)

c22.loc[r.index, 'c'] = c2.loc[r.index, 'c'].values.astype(int)


# or
#c22.loc[c22.val == c22.val.min()]


grid.to_csv(PATH + 'grid') # intermittent save

c22.to_csv(PATH + 'els')

# read back files
#grid = pd.read_csv('grid',index_col=0)


# Mask Land

#Define the coastline shapefile
shapefile = '/Users/brey/DATA/COASTLINES/naturalearth/coastline/ne_%sm_coastline' %{'l':110, 'i':50, 'h':10}['i']

#read it into a DataFrame
shp = gp.GeoDataFrame.from_file(shapefile+'.shp')

shp['length']=shp['geometry'][:].length # optional

shp = shp.sort_values(by='length', ascending=0) #optional
shp = shp.reset_index(drop=True)

#put all Lines in a list
ls=[]
for i in range(shp.shape[0]):
    il = shp.loc[i,'geometry']
    try:
        print( len(il))
        for k in range(len(list(il.geoms))):
               ls.append(list(il.geoms)[k])
    except:
        ls.append(il)

sall = geometry.MultiLineString(ls) #join them into a Multiline

c = ops.linemerge(sall) #merge parts if possible

print len(c)

# find the open LineStrings
op=[]
for i in range(len(c)):
    if not c[i].is_ring : op.append(i)

#compute the min/max latitude of the shapefile
latmin = 100.
latmax = -100.
for line in c:
    ylat = np.array([lat for (lon,lat) in line.coords[:]])
    latmin = min([latmin, ylat.min()])
    latmax = max([latmax, ylat.max()])


print latmin, latmax

points = grid.loc[:,['lon','lat']].values

# for all closed lines, select nodes inside
xi=[]
yi=[]
for i in range(len(c)):
    z = geometry.Polygon(c[i])
    path = mpltPath.Path(list(zip(z.boundary.xy[0],z.boundary.xy[1])))

    inside = path.contains_points(points)

#    if np.sum(inside) > 0:
    X = np.ma.masked_array(xp,mask=np.invert(inside)),
    Y = np.ma.masked_array(yp,mask=np.invert(inside))
    xi.append(X)
    yi.append(Y)

print len(xi)

#merge the masks 
gmask=np.ones(xi[0][0].shape, dtype=bool)
for i in range(len(xi)):
    gmask = np.logical_and(gmask,xi[i][0].mask)

#Mask also the south pole area
ylat = np.array([lat for (lon,lat) in c[0].coords[:]])
#points[points[:,1] < latmin ]
gmask[points[:,1] <= ylat[0]] = False

#dry points
gx = points[~gmask][:,0]
gy = points[~gmask][:,1]

# wet points
wx = points[gmask][:,0]
wy = points[gmask][:,1]


np.savez(PATH + 'gmask',gmask) # save mask


# identify the node numbers of land points
idx=[]
for i in range(len(gx)):
    idx.append(grid.loc[(grid['lon'] == gx[i]) & (grid['lat'] == gy[i])].index.values)

land = np.array([hf for hf in list(idx)]).flatten()


np.savez(PATH + 'land_index',land) # save 


# identify the node numbers of land points
iland=[]
for i in range(len(xi)):
    ibx = points[~xi[i][0].mask][:,0]# land nodes for the specific closed? shapefile
    iby = points[~xi[i][0].mask][:,1]
    print (i, len(ibx))
    bdx=[]
    for j in range(len(ibx)):
        bdx.append(grid.loc[(grid['lon'] == ibx[j]) & (grid['lat'] == iby[j])].index.values)
    iland.append(np.array([hf for hf in list(bdx)]).flatten())


bg=grid.copy()
be=c22.copy()


# ### function

bns=[]
for i in range(len(xi)):
    lnodes = bg.loc[iland[i]]
    els=[]
    if lnodes.size == 0: 
        bns.append([])
        continue
    for node in lnodes.index:
       # print('Node {}'.format(node))
        eln = be.loc[(be['a'] == node) | (be['b'] == node) | (be['c'] == node)] #find elemnets with these nodes
        els.append(eln.index.values)
    els = np.unique(np.hstack(els))
    b=[]
    for el in els:
        b.append(be.loc[el,['a','b','c']].values)
    b = np.unique(b)
    b = b.astype(int)
    hb=[]
    for bn in b:
        if bn not in lnodes.index.values : hb.append(bn)
    bns.append(hb)


np.savez(PATH + 'bns',bns) # save 


# Crop land 

gcrop=grid.copy()
ecrop=c22.copy()

#For the land nodes
lnodes = gcrop.loc[land]


els=[]
for node in lnodes.index:
   # print('Node {}'.format(node))
    eln = ecrop.loc[(ecrop['a'] == node) | (ecrop['b'] == node) | (ecrop['c'] == node)] #find elemnets with these nodes
    els.append(eln.index.values) # store the element index

    gcrop = gcrop.drop(node) # drop node

gcrop.reset_index(inplace=True) #keep the old index?

gcrop.index = gcrop.index + 1 # add one so that it starts from 1

dropels = list(np.hstack(els)) #there are the elements to be dropped - the ones which have dropped nodes

dropels = np.unique(dropels)


# Re - index original grid to exclude land elements

#for i in range(len(lnodes.index)):
#    node = lnodes.index[i] - i # everytime we crop the original index is reduced by one
#    eln = ecrop.loc[(ecrop['a'] == node) | (ecrop['b'] == node) | (ecrop['c'] == node)] #find elemnets with these nodes

#    ecrop.loc[ecrop['a']>=node,'a'] -= 1 # subtract one from all nodes > node  
#    ecrop.loc[ecrop['b']>=node,'b'] -= 1 
#    ecrop.loc[ecrop['c']>=node,'c'] -= 1 
        
#    ecrop = ecrop.drop(eln.index)
#    ecrop.reset_index(inplace=True, drop=True)

#vectorize the above???    
for i in range(len(lnodes.index)):
    node = lnodes.index[i] - i # everytime we crop the original index is reduced by one
    emask = (ecrop.a.values == node) | (ecrop.b.values == node) | (ecrop.c.values == node) #find elemnets with these nodes
    eln = ecrop[emask]

    amask = ecrop.a.values >= node
    ecrop.loc[amask,'a'] -= 1
    bmask = ecrop.b.values >= node
    ecrop.loc[bmask,'b'] -= 1
    cmask = ecrop.c.values >= node
    ecrop.loc[cmask,'c'] -= 1

    ecrop = ecrop.drop(eln.index)
    ecrop.reset_index(inplace=True, drop=True)

ecrop.index = ecrop.index + 1


### find hanged nodes (leftover nodes without elements)
tri = ecrop.loc[:,['a','b','c']].values - 1 # for the plot below (index.min => 0)
q = np.unique(tri.flatten()+1) # all the unique nodes in elements
dq = list(set(range(1,gcrop.shape[0])) - set(q)) # the ones that are in gcrop but not in elems
dq.sort()


#drop them from gcrop 
for i in range(len(dq)):
    node = dq[i] - i # everytime we crop the original index is reduced by one
    eln = ecrop.loc[(ecrop['a'] == node) | (ecrop['b'] == node) | (ecrop['c'] == node)] #find elemnets with these nodes

    ecrop.loc[ecrop['a']>=node,'a'] -= 1 # subtract one from all nodes > node  
    ecrop.loc[ecrop['b']>=node,'b'] -= 1 
    ecrop.loc[ecrop['c']>=node,'c'] -= 1 
        
    ecrop = ecrop.drop(eln.index)
    ecrop.reset_index(inplace=True, drop=True)

ecrop.index = ecrop.index + 1

gcrop = gcrop.drop(dq)

gcrop.reset_index(inplace=True, drop=True)

gcrop.index = gcrop.index + 1

tri = ecrop.loc[:,['a','b','c']].values - 1 # for the plot below (index.min => 0)

#rename index column WE NEED IT FOR THE BOUNDARIES BELOW
gcrop.columns = ['bindex', 'x', 'y', 'z', 'q', 'lon', 'lat']


# SAVE to file for import
gcrop.to_csv(PATH + 'gcrop')
ecrop.to_csv(PATH + 'ecrop')

#gcrop = pd.read_csv('gcrop',index_col=0)
#ecrop = pd.read_csv('ecrop',index_col=0)

### Read Bathymetry

fbath = '/Users/brey/GitHub/jigsaw-geo-matlab/jigsaw/geo/topo.msh'
coord, nc = pd.read_csv(fbath,header=2,skiprows=1,nrows=0,delimiter='=')
coord, nc
nc0, nc1 = nc.split(';')
nlon = int(nc1)
#read lon for all boundary nodes
xb = pd.read_csv(fbath,skiprows=4,header=None,nrows=nlon,names=['lon'])
coord, nc = pd.read_csv(fbath,header=4,skiprows=nlon,nrows=0,delimiter='=')
coord, nc
nc0, nc2 = nc.split(';')
nlat = int(nc2)
#read lat,depth for boundary nodes
yb = pd.read_csv(fbath,skiprows=5 + nlon,header=None,nrows=nlat,names=['lat'])
value, nc = pd.read_csv(fbath,header=4,skiprows=nlon+nlat+1,nrows=0,delimiter='=')
value, nc
val, val0 = nc.split(';')
val = int(val)
#read 
bath = pd.read_csv(fbath,skiprows=6 + nlon + nlat,header=None,nrows=val,names=['value'])
vb = bath.values.reshape(nlon,nlat)
zlev = vb.T


### resample on grid points
X, Y = np.meshgrid(xb.lon.values, yb.lat.values)

#mask the positive values
wet = np.ma.masked_array(zlev,zlev>0)
# wet.fill_value = 0.
mx = np.ma.masked_array(X,wet.mask) 
my = np.ma.masked_array(Y,wet.mask)


orig = pyresample.geometry.SwathDefinition(lons=mx,lats=my) # original points
targ = pyresample.geometry.SwathDefinition(lons=gcrop.lon.values,lats=gcrop.lat.values) # target grid
       
ibath = pyresample.kd_tree.resample_nearest(orig,wet,targ,radius_of_influence=50000,fill_value=0)

gcrop['dep'] = -ibath # make bathymetry positive 

nn = gcrop.shape[0]
ne = ecrop.shape[0]


gcrop.to_csv(PATH + 'gcrop_b')


### save to a gr3 file

### Compute internal land boundaries

#### All ...

lbound = [x for x in bns if len(x) > 0]

bound=pd.DataFrame({})
for i in range(len(lbound)):
        bound = pd.concat([bound , gcrop.loc[gcrop['bindex'].isin(lbound[i])]])
     
bound.to_csv(PATH + 'bound')


# ### Combine points

## Find the connections with other lists of boundary points

ibs=[]
for i in range(len(lbound)):
        ibs.append(list(gcrop.loc[gcrop['bindex'].isin(lbound[i])].index.values))

idb=ibs[0]
for i in range(1,len(ibs)):
    ibels = bels.loc[(ecrop.a.isin(ibs[i])) | (bels.b.isin(ibs[i])) | (bels.c.isin(ibs[i]))] # all other nodes that share elements with the lbs[0] nodes 
    if ibels.shape[0] > 0 :
#        qq = ibels.loc[(ibels.a.isin(ibs[0])) | (ibels.a.isin(ibs[0])) | (ibels.a.isin(ibs[0]))] # all other nodes that share nodes with the lbs[0] nodes
#        if qq.shape[0] > 0:
            idb = idb + ibs[i]

#### output

folder = '/Users/brey/GLOBAL/SCHISM_GLOBAL_2/'

g3file = folder + 'hgrid.gr3'

with open(g3file,'w') as f:
    f.write('\t uniform.gr3\n')
    f.write('\t {} {}\n'.format(ne,nn))
    
gcrop.to_csv(g3file,index=True, sep='\t', header=None,mode='a', float_format='%.10f', columns=['lon','lat','dep'])

ecrop.to_csv(g3file,index=True, sep='\t', header=None, mode='a', columns=['q','a','b','c'])

with open(g3file, 'a') as f:
    f.write('{} = Number of open boundaries\n'.format(0))
    f.write('{} = Total number of open boundary nodes\n'.format(0))

with open(g3file, 'a') as f:
    f.write('{} = Number of land boundaries\n'.format(len(ibs)))
    f.write('{} = Total number of land boundary nodes\n'.format(bound.shape[0]))

    total = 0 
    for i in range(len(ibs)):
        f.write('{} 0 = Number of nodes for land boundary {}\n'.format(len(ibs[i]), i+1))
        total += len(ibs[i])
        for item in ibs[i]:
            f.write("%s\n" % item)

llfile=folder+'hgrid.ll'

copyfile(g3file, llfile)

manfile=folder+'manning.gr3'

with open(manfile,'w') as f:
    f.write('\t 0 \n')
    f.write('\t {} {}\n'.format(ne,nn))


gcrop['man']=.12
gcrop.head()

gcrop.to_csv(manfile,index=True, sep='\t', header=None,mode='a', float_format='%.10f',columns=['lon','lat','man'] )

ecrop.to_csv(manfile,index=True, sep='\t', header=None, mode='a', columns=['q','a','b','c'])

