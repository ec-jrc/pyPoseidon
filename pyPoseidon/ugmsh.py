"""
gmsh module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import pandas as pd
import numpy as np
import geopandas as gp
import xarray as xr
import os
from tqdm import tqdm
import sys
import gmsh

from pyPoseidon.utils.tag import *
from pyPoseidon.utils.spline import *
from pyPoseidon.utils.stereo import to_lat_lon

import logging
logger = logging.getLogger('pyPoseidon')


def read_gmsh(mesh,**kwargs):

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()

    gmsh.open(mesh)

    logger.info('Analyze grid')

    nodeTags, coord, parametricCoord = gmsh.model.mesh.getNodes()

    nodes = pd.DataFrame(coord.reshape(-1,3),columns=['x','y','z'])

    elementTags2, nodeTags2 = gmsh.model.mesh.getElementsByType(2)

    elems = nodeTags2.reshape(-1,3)

    tria = pd.DataFrame(elems-1,columns=['a','b','c'])

    # boundaries

    bounds=[]

    bgs=gmsh.model.getPhysicalGroups()

    bgs=pd.DataFrame(bgs[:-1],columns=['dim','tag'])

    #open boundaries
    logger.info('open boundaries')

    obs=bgs.loc[bgs.tag<1000]

    for row in obs.itertuples(index=True, name='Pandas'):
        onodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"),tag=getattr(row, "tag"))

        db = pd.DataFrame({'node':onodes-1})
        db['type']=np.nan
        db['id']=getattr(row, "Index")+1

        bounds.append(db)


    #land boundaries type
    logger.info('land boundaries')

    lbs=bgs.loc[(bgs.tag>1000)&(bgs.tag<2000)]
    lbs.reset_index(inplace=True,drop=True)

    for row in lbs.itertuples(index=True, name='Pandas'):
        lnodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"),tag=getattr(row, "tag"))

        db = pd.DataFrame({'node':lnodes-1})
        db['type']=0
        db['id']=-(getattr(row, "Index")+1)

        bounds.append(db)

    if lbs.empty: #Store max index from land boundaries
        itag=0
    else:
        itag=lbs.tag.max()-1000


    #islands
    logger.info('islands')


    ibs=bgs.loc[bgs.tag>2000]
    ibs.reset_index(inplace=True,drop=True)
    ibs.index=ibs.index + itag #set index


    for row in ibs.itertuples(index=True, name='Pandas'):

        inodes, xyz = gmsh.model.mesh.getNodesForPhysicalGroup(dim=getattr(row, "dim"),tag=getattr(row, "tag"))
        db = pd.DataFrame({'node':inodes-1})
        db['type']=-1
        db['id']=-(getattr(row, "Index")+1)

        bounds.append(db)

    if bounds != []:
        bnodes=pd.concat(bounds).reset_index(drop=True)

        bnodes.index.name='bnodes'

        bnodes = bnodes.drop_duplicates('node')

        bnodes['id']=bnodes.id.astype(int)

    else:
        bnodes=pd.DataFrame({})


    #check if global and reproject
    sproj = kwargs.get('gglobal', False)
    if sproj : # convert to lat/lon
        xd,yd = to_lat_lon(nodes.x,nodes.y)
        nodes['x']=xd
        nodes['y']=yd

    grid = pd.DataFrame({'lon':nodes.x, 'lat':nodes.y})

    tri3 = tria.values

    logger.info('Finalize Dataset')

    ## make dataset
    els = xr.DataArray(
          tri3,
          dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'], name='SCHISM_hgrid_face_nodes'
          )

    nod = grid.loc[:,['lon','lat']].to_xarray().rename({'index':'nSCHISM_hgrid_node','lon':'SCHISM_hgrid_node_x','lat':'SCHISM_hgrid_node_y'})
    nod = nod.drop_vars('nSCHISM_hgrid_node')

    dep = xr.Dataset({'depth': (['nSCHISM_hgrid_node'], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    gr = xr.merge([nod,els,dep,bnodes.to_xarray()])

    gmsh.finalize()

    return gr


def gmsh_(**kwargs):

    logger.info('Creating grid with GMSH\n')

    geometry = kwargs.get('geometry', None)

    if isinstance(geometry,dict):

        rpath = kwargs.get('rpath', '.')

        if not os.path.exists(rpath):
                os.makedirs(rpath)

        gpath=os.path.join(rpath,'gmsh')
        if not os.path.exists(gpath):
                os.makedirs(gpath)

        df , bmindx = tag_(**kwargs)

        make_gmsh(df, **kwargs)

        gr = read_gmsh(rpath+'/gmsh/mymesh.msh')

    return gr



def gset(df,**kwargs):

    logger.info('interpolate coastal points')

    lc = kwargs.get('lc', .5)

    df['lc'] = lc
    df = df.apply(pd.to_numeric)

    #Resample to equidistant points

    conts = np.unique(df.index[df.tag<0].get_level_values(0))
    conts = [x for x in conts if x not in ['line0']]#except the outer LineString

    ibs=len(conts)

    ndfsfs={}
    for ic in tqdm(range(ibs)):
        contour = conts[ic]
        curve = df.loc[contour,['lon','lat']]
        curve = pd.concat([curve,curve.loc[0:0]]).reset_index(drop=True)
        di = spline(curve,ds=.01,method='slinear')
        di['z']=df.loc[contour].z.values[0]
        di['tag']=df.loc[contour].tag.values[0].astype(int)
        di['lc']=df.loc[contour].lc.values[0]
        ndfsfs.update({contour:di.drop_duplicates(['lon','lat'])})

    df_ = pd.concat(ndfsfs, axis=0)
    df_['z'] = df_.z.values.astype(int)
    df_['tag'] = df_.tag.values.astype(int)

    #Line0

    logger.info('set outermost boundary')

    df0=df.loc['line0']

    mtag=df0.tag.min()
    mtag=mtag.astype(int)

    nd0={}
    for ic in tqdm(range(mtag,0)):
        contour = df0.tag==ic
        curve = df0.loc[contour,['lon','lat']].reset_index(drop=True)
    #    curve = pd.concat([curve,curve.loc[0:0]]).reset_index(drop=True)
        di = spline(curve,ds=.01,method='slinear')
        di['z']=df0.loc[contour].z.values[0]
        di['tag']=df0.loc[contour].tag.values[0].astype(int)
        di['lc']=df0.loc[contour].lc.values[0]
        nd0.update({ic:di.drop_duplicates(['lon','lat'])})

    # Join Line0
    df0_=df0.copy()
    for l in range(mtag,0):
    #    print(l)
        idx=df0_.loc[df0_.tag==l].index
        df0_=pd.concat([df0_.iloc[:idx[0]],nd0[l],df0_.iloc[idx[-1]+1:]])
        df0_.reset_index(drop=True, inplace=True)

    df0_=pd.concat({'line0': df0_})

    #join all
    ddf = pd.concat([df0_,df_])

    ddf['z'] = ddf.z.values.astype(int)
    ddf['tag'] = ddf.tag.values.astype(int)

    #check orientation
    r0=ddf.loc['line0']

    if not shapely.geometry.LinearRing(r0[['lon','lat']].values).is_ccw:

        rf0=ddf.loc['line0'].iloc[::-1].reset_index(drop=True)
        ddf.loc['line0']=rf0.values

    ddf = ddf.apply(pd.to_numeric)

    return ddf


def make_gmsh(df, **kwargs):

    logger.info('create grid')

    model = gmsh.model
    factory = model.geo

    gmsh.initialize()
    model.add("schism")

#    gmsh.option.setNumber("General.Terminal", 1)

    interpolate = kwargs.get('interpolate', False)
    if interpolate:
        ddf=gset(df,**kwargs)
    else:
        ddf=df
    lc = kwargs.get('lc', .5)

    ddf['lc'] = lc
    ddf = ddf.apply(pd.to_numeric)


    # save boundary configuration for Line0
    rb0=ddf.loc['line0'].copy()

    if not shapely.geometry.LinearRing(rb0[['lon','lat']].values).is_ccw:# check for clockwise orientation
        rb0=ddf.loc['line0'].iloc[::-1].reset_index(drop=True)

    rb0.index=rb0.index+1 # fix index
    rb0['bounds']=[[i,i+1] for i in rb0.index]
    rb0['bounds']=rb0.bounds.values.tolist()[:-1]+[[rb0.index[-1],1]] # fix last one

    #store blines
    blines={}

    for tag_ in rb0.tag.unique():

        ibs=rb0.loc[rb0.tag==tag_].index.values
        #ibs

        lbs=rb0.loc[rb0.tag==tag_].bounds.values.tolist()
        #lbs

        ai=np.unique(np.concatenate(lbs))
        #ai

        itags=[i for i in ai if i in ibs]
        #itags

        if tag_ > 0:
            items = set(itags)

            imask=[set(x).issubset(items) for x in rb0.loc[rb0.tag==tag_].bounds]
            #imask

            bi=rb0.loc[rb0.tag==tag_][imask].index.values.tolist()


        else:

            bi=rb0.loc[rb0.tag==tag_].index.values.tolist()

        blines.update({tag_:bi})

    al = [j for i in list(blines.values()) for j in i]
    lover=[x for x in rb0.index if x not in al]

    for i,v in rb0.loc[lover].iterrows():
        nns=rb0.loc[v[5],['tag']].values
        itag=[x for x in nns if x < 0 ][0]
        blines.update({itag[0]:blines[itag[0]]+[i]})


    land_lines = { your_key: blines[your_key] for your_key in [x for x in blines.keys() if x < 0] }
    open_lines = { your_key: blines[your_key] for your_key in [x for x in blines.keys() if x > 0] }

    logger.info('Define geometry')

    loops=[]
    islands=[]
    all_lines=[]

    ltag=1

    for row in rb0.itertuples(index=True, name='Pandas'):
        factory.addPoint(getattr(row, "lon"),getattr(row, "lat"),getattr(row, "z"),getattr(row, "lc"),getattr(row, "Index"))
    for row in rb0.itertuples(index=True, name='Pandas'):
        factory.addLine(getattr(row, "bounds")[0],getattr(row, "bounds")[1],getattr(row, "Index"))

    lines=rb0.index.values
    all_lines.append(lines)

    tag=rb0.index.values[-1]

    factory.addCurveLoop(lines, tag=ltag)
    #print(loop)
    loops.append(ltag)
    all_lines.append(lines)

    tag += 1
    ltag += 1

    for contour in tqdm(ddf.index.levels[0][1:]):
        rb=ddf.loc[contour].copy()
        if not shapely.geometry.LinearRing(rb[['lon','lat']].values).is_ccw:# check for clockwise orientation
            rb=ddf.loc[contour].iloc[::-1].reset_index(drop=True)

        rb.index=rb.index+tag
        rb['bounds']=[[i,i+1] for i in rb.index]
        rb['bounds']=rb.bounds.values.tolist()[:-1]+[[rb.index[-1],rb.index[0]]] # fix last one


        for row in rb.itertuples(index=True, name='Pandas'):
            factory.addPoint(getattr(row, "lon"),getattr(row, "lat"),getattr(row, "z"),getattr(row, "lc"),getattr(row, "Index"))
        for row in rb.itertuples(index=True, name='Pandas'):
            factory.addLine(getattr(row, "bounds")[0],getattr(row, "bounds")[1],getattr(row, "Index"))

        lines=rb.index.values
        all_lines.append(lines)

        tag=rb.index.values[-1]+1

        factory.addCurveLoop(lines,tag=ltag)
    #    print(tag)
        loops.append(ltag)

        islands.append(lines)
        all_lines.append(lines)


        tag += 1
        ltag += 1

    factory.addPlaneSurface(loops)
    logger.info('synchronize')
    factory.synchronize()

    ## Group open boundaries lines
    for key,values in open_lines.items():
        gmsh.model.addPhysicalGroup(1, values,1000-int(key))

    ## Group land boundaries lines
    for key,values in land_lines.items():
        gmsh.model.addPhysicalGroup(1, values,1000-int(key))

    ntag=1
    for k in tqdm(range(len(islands))):
        gmsh.model.addPhysicalGroup(1, islands[k], 2000+ntag)
        ntag += 1

    ps = gmsh.model.addPhysicalGroup(2, [1])
    gmsh.model.setPhysicalName(2, ps, "MyMesh")

    flat_list = [item for sublist in all_lines for item in sublist]
    ols=[j for i in list(open_lines.values()) for j in i]
    lists = [x for x in flat_list if x not in ols]

    model.mesh.field.add("Distance", 1)
    model.mesh.field.setNumbers(1, "CurvesList", lists)

    model.mesh.field.add("Threshold", 2);
    model.mesh.field.setNumber(2, "InField", 1);
    model.mesh.field.setNumber(2, "SizeMin", .01);
    model.mesh.field.setNumber(2, "SizeMax", .1);
    model.mesh.field.setNumber(2, "DistMin", .01);
    model.mesh.field.setNumber(2, "DistMax", .1);
#    model.mesh.field.setNumber(2, "StopAtDistMax", 1);

    # Merge a post-processing view containing the target anisotropic mesh sizes
#    gmsh.merge('iceland.pos')

#    model.mesh.field.add("PostView", 3)
#    model.mesh.field.setNumber(3, "ViewIndex", 0)

#    model.mesh.field.add("Min", 4)
#    model.mesh.field.setNumbers(4, "FieldsList", [2,3])

    model.mesh.field.setAsBackgroundMesh(2)

    gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
    gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0)
    gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0)

    logger.info('execute')

    gmsh.model.mesh.generate(2)

    # ... and save it to disk
    rpath = kwargs.get('rpath', '.')

    logger.info('save mesh')
#    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.write(rpath + '/gmsh/mymesh.msh')

#    gmsh.write('mymesh.vtk')

    gmsh.finalize()
