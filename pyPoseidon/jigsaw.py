"""
Jigsaw module

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
import shapely
import subprocess
import sys

from pyPoseidon.utils.stereo import to_lat_lon, to_stereo
from pyPoseidon.utils.sort import *
import pyPoseidon.dem as pdem
from pyPoseidon.utils.hfun import *
import logging        
        
logger = logging.getLogger('pyPoseidon')




def geo(df, path='.', tag='jigsaw'):
    
    fgeo = path + tag+'-geo.msh'
    # write header
    with open(fgeo,'w') as f:
        f.write('#{}; created by pyPoseidon\n'.format(tag +'-geo.msh'))
        f.write('MSHID=2;EUCLIDEAN-MESH\n')
        f.write('NDIMS=2\n')
        f.write('POINT={}\n'.format(df.shape[0]))
    
    #write lines
    with open(fgeo, 'a') as f:
        for line in df.index.levels[0]:
            df.loc[line].to_csv(f, index=False, header=0, columns=['lon','lat','z'],sep=';')


    edges = pd.DataFrame([]) # initiate

    # create edges
    for line in df.index.levels[0]:
        i0 = edges.shape[0]
        ie = df.loc[line].shape[0] + edges.shape[0]
        dline = df.loc[line].copy()
        dline.index = range(i0,ie)
        dline.loc[:,'ie'] = dline.index.values + 1
        dline.loc[dline.index[-1],'ie']=i0
        dout = dline.reset_index().loc[:,['index','ie','tag']]

        dout['tag1'] = dline.loc[dline.ie.values,'tag'].values.astype(int)
        dout['que'] = np.where(((dout['tag'] != dout['tag1']) & (dout['tag'] > 0)) , dout['tag1'], dout['tag'])    
        dout = dout.reset_index().loc[:,['index','ie','que']]
        edges = edges.append(dout)        
        
    # write header
    with open(fgeo,'a') as f:
        f.write('EDGE2={}\n'.format(edges.shape[0]))
    
    with open(fgeo, 'a') as f:
        edges.to_csv(f, index=False, header=0, sep=';')
        


def read_msh(fmsh):
    
    grid = pd.read_csv(fmsh, header=0, names=['data'], index_col=None, low_memory=False)
    npoints = int(grid.loc[2].str.split('=')[0][1])

    nodes = pd.DataFrame(grid.loc[3: 3 + npoints - 1,'data'].str.split(';').values.tolist(),columns=['lon','lat','z'])

    ie = grid[grid.data.str.contains('EDGE')].index.values[0]
    nedges = int(grid.loc[ie].str.split('=')[0][1])
    edges = pd.DataFrame(grid.loc[ie + 1 :ie  + nedges ,'data'].str.split(';').values.tolist(),columns=['e1','e2','e3'])

    i3 = grid[grid.data.str.contains('TRIA')].index.values[0]
    ntria = int(grid.loc[i3].str.split('=')[0][1])
    tria = pd.DataFrame(grid.loc[i3 + 1 : i3 + ntria + 1 ,'data'].str.split(';').values.tolist(),columns=['a','b','c','d'])

    return [nodes,edges,tria]
    
    
            

def jigsaw(**kwargs):
     
    
    logger.info('Creating grid with JIGSAW\n')
       
    geometry = kwargs.get('geometry', None)
    
    if isinstance(geometry,dict): 
             
        df , bmindx = jdefault(**kwargs)
        
        gr = jigsaw_(df, bmindx, **kwargs)       
        
    elif isinstance(geometry,str):
        
        if geometry=='global':
                 
            hfun0 = hfun_(kwargs.get('coastlines',None),kwargs.get('res',.1))
            
            kwargs.update({'hfun':hfun0}) 
            
            df = sgl(**kwargs)
            
            bmindx = df.tag.min()
            
            gr = jigsaw_(df, bmindx, **kwargs)   
            
            
            convert = kwargs.get('to_lat_lon',True)
            
            if convert :
                # convert to lat/lon
                u, v = gr.SCHISM_hgrid_node_x.values, gr.SCHISM_hgrid_node_y.values 
            
                rlon, rlat = to_lat_lon(u,v)
            
                gr['SCHISM_hgrid_node_x'].values = rlon
                gr['SCHISM_hgrid_node_y'].values = rlat
                       
        else:
        
            df = jcustom(**kwargs)
        
            bmindx = df.tag.min()
        
            gr = jigsaw_(df, bmindx, **kwargs)
          
    
    return gr


def sgl(**kwargs):
    
    try:
        geo = gp.GeoDataFrame.from_file(kwargs.get('coastlines',None))
    except:
        logger.warning('coastlines is not a file, trying with geopandas Dataset')
        try:
            geo = kwargs.get('coastlines',None)
        except:    
            logger.error('coastlines argument not valid ')
            sys.exit(1)
    
        # Manage coastlines
        logger.info('preparing coastlines')
        #Kaspian Sea, if present
        kasp = ((geo.bounds.miny > 36.) & (geo.bounds.maxy < 48.) & (geo.bounds.maxx < 55.) & (geo.bounds.minx > 45.))
        geo = geo.drop(geo.loc[kasp].index)
    
        #ANTARTICA
        anta_mask = (geo.bounds.miny < -89.) # indentify antartica
        anta = geo.loc[anta_mask]
        indx = anta.index # keep index
    
        anta = pd.DataFrame(anta.boundary.values[0].coords[:], columns=['lon','lat']) # convert boundary values to pandas
        d1 = anta.where(anta.lon==anta.lon.max()).dropna().index[1:] # get artificial boundaries as -180/180
        d2 = anta.where(anta.lon==anta.lon.min()).dropna().index[1:]
        anta = anta.drop(d1).drop(d2) # drop the points
        d3 = anta.where(anta.lat==anta.lat.min()).dropna().index # drop lat=-90 line
        anta = anta.drop(d3)
        an = gp.GeoDataFrame({'geometry' : [shapely.geometry.LineString(anta.values)], 'length':shapely.geometry.LineString(anta.values).length }, index=indx) # put together a LineString
        geo.loc[indx] = shapely.geometry.LineString(anta.values) # put it back to geo
    
        # International Meridian
        m1 = geo[geo.bounds.minx == -180.].index
        m2 = geo[geo.bounds.maxx == 180.].index
        mm = np.concatenate((m1, m2)) # join them
        mm = [j for j in mm if j!=indx] # subtract antartica
    
        # convert to u,v (stereographic coordinates)
        for idx, poly in geo.iterrows():
            geo.loc[idx,'geometry'] = shapely.ops.transform(lambda x,y,z=None: to_stereo(x,y), poly.geometry)
    
        w = geo.drop(indx) # get all polygons
        ww = w.loc[mm] # join the split polygons
        gw = gp.GeoDataFrame(geometry = list(ww.buffer(.0001).unary_union)) # merge the polygons that are split (around -180/180)
    
        w = w.drop(mm)
        # Check antartica LineString
        if not geo.iloc[indx].geometry.values[0].is_ring:
            ca = gp.GeoDataFrame(geometry = [shapely.geometry.LinearRing(geo.loc[indx].geometry.values[0])], index=indx)
            ca['geometry'] = shapely.geometry.LineString(ca.geometry.values[0])
        else:
            ca = geo.loc[indx]
    
        # PUT ALL TOGETHER    
        geo = pd.concat([w, gw, ca], ignore_index=True).reset_index(drop=True)        
    
    logger.info('storing boundaries')
    
    geo['tag'] = - (geo.index + 1)
    
    idx=0
    dic={}
    for i, line  in geo.iloc[:-1].iterrows():
        lon=[]
        lat=[]
        try:
            for x,y in line.geometry.boundary.coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({'line{}'.format(idx):{'lon':lon,'lat':lat,'tag':line.tag}})
            idx += 1
        except:
            for x,y in line.geometry.boundary[0].coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({'line{}'.format(idx):{'lon':lon,'lat':lat,'tag':line.tag}})
            idx += 1
            
    for i, line  in geo.iloc[-1:].iterrows():
        lon=[]
        lat=[]
        for x,y in line.geometry.coords[:]:
            lon.append(x)
            lat.append(y)
            dic.update({'line{}'.format(idx):{'lon':lon,'lat':lat,'tag':line.tag}})

    dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()}

    df = pd.concat(dict_of_df, axis=0)

    df['z']=0
    df = df.drop_duplicates() # drop the repeat value on closed boundaries

    return df


def jcustom(**kwargs):


    geometry = kwargs.get('geometry', None)

    try:
        geo = gp.GeoDataFrame.from_file(geometry)
    except:
        logger.error('geometry argument not a valid file')
        sys.exit(1)


    idx=0
    dic={}
    for idx, line  in geo.iterrows():
        lon=[]
        lat=[]
        for x,y in line.geometry.coords[:]:
            lon.append(x)
            lat.append(y)
        dic.update({'line{}'.format(idx):{'lon':lon,'lat':lat,'tag':line.tag}})
        idx += 1

    dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()}
    
    df = pd.concat(dict_of_df, axis=0)
    
    df['z']=0
    df = df.drop_duplicates() # drop the repeat value on closed boundaries
    
    out_b = []
    for line in df.index.levels[0]:
        out_b.append(shapely.geometry.LineString(df.loc[line,['lon','lat']].values))

    merged = shapely.ops.linemerge(out_b)
    merged = pd.DataFrame(merged.coords[:], columns=['lon','lat'])
    merged = merged.drop_duplicates()
    match = df.drop_duplicates(['lon','lat']).droplevel(0)
    match = match.reset_index(drop=True)

    df1 = merged.sort_values(['lon', 'lat'])
    df2 = match.sort_values(['lon', 'lat'])
    df2.index = df1.index
    final = df2.sort_index()
    final = pd.concat([final],keys=['line0'])

    
    return final


def jdefault(**kwargs):
    
    world = kwargs.get('coastlines',None)
    
    if world is None :
        logger.error('coastlines not given')
        sys.exit(1)
    
    world = world.explode()
    
    geometry = kwargs.get('geometry', None)
    
    try:
        lon_min = geometry['lon_min']
        lon_max = geometry['lon_max']
        lat_min = geometry['lat_min']
        lat_max = geometry['lat_max']
    except:
        logger.error('geometry not set properly')
        sys.exit(1)
    
    
    #create a polygon of the lat/lon window
    grp=shapely.geometry.Polygon([(lon_min,lat_min),(lon_min,lat_max),(lon_max,lat_max),(lon_max,lat_min)])

    #create a LineString of the grid
    grl=shapely.geometry.LineString([(lon_min,lat_min),(lon_min,lat_max),(lon_max,lat_max),(lon_max,lat_min),(lon_min,lat_min)])
    
    # check -180/180 trespass
    if np.mean([lon_min,lon_max]) < 0 and lon_min < -180. :
        flag = -1
    elif np.mean([lon_min,lon_max]) > 0 and lon_max > 180. :
        flag = 1
    else:
        flag = 0

    #adjust abd mask based on lat/lon window
    if flag == 1 :
        block1 = world.cx[lon_min:180,lat_min:lat_max].copy()
        block2 = world.cx[-180:(lon_max-360.),lat_min:lat_max].copy()
    
        for idx, poly in block2.iterrows():
            block2.loc[idx,'geometry'] = shapely.ops.transform(lambda x,y,z=None: (x+360.,y), poly.geometry)
    
        block = block1.append(block2)
    
    elif flag == -1 :
    
        block1 = world.cx[lon_min + 360 : 180,lat_min:lat_max].copy()
        block2 = world.cx[-180:lon_max,lat_min:lat_max].copy()
    
        for idx, poly in block1.iterrows():
            block1.loc[idx,'geometry'] = shapely.ops.transform(lambda x,y,z=None: (x - 360.,y), poly.geometry)
    
        block = block1.append(block2)

    else:
        block = world.cx[lon_min:lon_max,lat_min:lat_max]

        
    g = block.buffer(.001).unary_union.symmetric_difference(grp) # get the dif from the world

    try: # make geoDataFrame 
        t = gp.GeoDataFrame({'geometry':g})
    except:
        t = gp.GeoDataFrame({'geometry':[g]})    
    t['length']=t['geometry'][:].length #get length
    t = t.sort_values(by='length', ascending=0) #sort
    t = t.reset_index(drop=True)
    
    t['in'] = gp.GeoDataFrame(geometry=[grp] * t.shape[0]).contains(t) # find the largest of boundaries
    idx = np.where(t['in']==True)[0][0] # first(largest) boundary within lat/lon
    b = t.iloc[idx].geometry #get the largest 
    
    # SETUP JIGSAW
    
    dic={}
    try:
        for l in range(len(b.boundary)):
            lon=[]
            lat=[]
            for x,y in b.boundary[l].coords[:]: 
                lon.append(x)
                lat.append(y)
            dic.update({'line{}'.format(l):{'lon':lon,'lat':lat}})
    except:
            lon=[]
            lat=[]
            for x,y in b.boundary.coords[:]: 
                lon.append(x)
                lat.append(y)
            dic.update({'line{}'.format(0):{'lon':lon,'lat':lat}})

    
    dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()}
    df = pd.concat(dict_of_df, axis=0)
    df['z']=0
    df = df.drop_duplicates() # drop the repeat value on closed boundaries
    
    
    #open (water) boundaries
    
    try:
        water = b.boundary[0] - (b.boundary[0] - grl)
    except:
        water = b.boundary - (b.boundary - grl)

    try:
        cwater = shapely.ops.linemerge(water)
    except:
        cwater= water

    mindx = 0
    #get all lines in a pandas DataFrame
    if cwater.type == 'LineString':
            lon=[]
            lat=[]
            for x,y in cwater.coords[:]:
                lon.append(x)
                lat.append(y)
            dic= {'line{}'.format(mindx):{'lon':lon,'lat':lat, 'z':0 ,'tag':mindx + 1}}
            mindx += 1
        
    elif cwater.type == 'MultiLineString' :
        dic={}
        for l in range(len(cwater)):
            lon=[]
            lat=[]
            for x,y in cwater[l].coords[:]:
                lon.append(x)
                lat.append(y)
            dic.update({'line{}'.format(mindx):{'lon':lon,'lat':lat, 'z':0 ,'tag':mindx + 1}})
            mindx += 1
    
    dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()}
    df_water = pd.concat(dict_of_df, axis=0)

       
    # land boundaries!!
    try:
        land = b.boundary[0] - grl
    except:
        land = b.boundary - grl 

    try:
        cland = shapely.ops.linemerge(land)
    except:
        cland = land
        
    mindx = 0 
    
    #get all lines in a pandas DataFrame
    dic_land = {}
    
    if cland.type == 'LineString':
            lon=[]
            lat=[]
            for x,y in cland.coords[:]:
                lon.append(x)
                lat.append(y)
            dic_land= {'line{}'.format(mindx):{'lon':lon,'lat':lat, 'z':0 ,'tag':mindx - 1}}
            mindx -= 1
        
    elif cland.type == 'MultiLineString' :
        dic_land={}
        for l in range(len(cland)):
            lon=[]
            lat=[]
            for x,y in cland[l].coords[:]:
                lon.append(x)
                lat.append(y)
            dic_land.update({'line{}'.format(mindx):{'lon':lon,'lat':lat, 'z':0 ,'tag':mindx - 1}})
            mindx -= 1
    
    dict_of_df = {k: pd.DataFrame(v) for k,v in dic_land.items()}

    try:       
        df_land = pd.concat(dict_of_df, axis=0)
    except:
        df_land = pd.DataFrame({})
    
    ddf = pd.concat([df_water,df_land])
    
    #Sort outer boundary
    
    out_b = []
    for line in ddf.index.levels[0]:
        out_b.append(shapely.geometry.LineString(ddf.loc[line,['lon','lat']].values))

    merged = shapely.ops.linemerge(out_b)
    merged = pd.DataFrame(merged.coords[:], columns=['lon','lat'])
    merged = merged.drop_duplicates()
    match = ddf.drop_duplicates(['lon','lat']).droplevel(0)
    match = match.reset_index(drop=True)

    df1 = merged.sort_values(['lon', 'lat'])
    df2 = match.sort_values(['lon', 'lat'])
    df2.index = df1.index
    final = df2.sort_index()
    final = pd.concat([final],keys=['line0'])
    
    
    bmindx = mindx
    
    #merge with islands    
    ndf = df.drop('line0')
    ndf['tag']=''
    for line in ndf.index.levels[0][1:]:
        ndf.loc[line,'tag'] = mindx - 1
        mindx -= 1
    ndf['tag'] = ndf.tag.astype(int)
    df = pd.concat([final,ndf])
    
    return df, bmindx      
    
def jigsaw_(df, bmindx, **kwargs):    
    
    logger.info('Creating JIGSAW files\n')
    
    tag = kwargs.get('tag', 'jigsaw')
    rpath = kwargs.get('rpath', '.')
    
    
    if not os.path.exists(rpath):
            os.makedirs(rpath)
    
    path = rpath+'/jigsaw/'
    if not os.path.exists(path):
        os.makedirs(path)
    
    
    geo(df,path=path,tag=tag)
          
    hfun = kwargs.get('hfun', None)
    
    if hfun is not None:

        if isinstance(hfun,str):
            dh = xr.open_dataset(hfun)
            to_hfun_grid(dh,path + tag+'-hfun.msh')   # write hfun file  
        else:
            to_hfun_mesh(hfun,path + tag+'-hfun.msh') 
    
    # write jig file
    fjig = path + '/' + tag+'.jig'
    
    with open(fjig,'w') as f:
        f.write('GEOM_FILE ={}\n'.format(tag+'-geo.msh'))
        f.write('MESH_FILE ={}\n'.format(tag+'.msh'))
        if hfun : f.write('HFUN_FILE ={}\n'.format(tag+'-hfun.msh'))
        f.write('HFUN_SCAL = ABSOLUTE\n')
        f.write('HFUN_HMAX = Inf\n')
        f.write('HFUN_HMIN = 0.0\n')
        f.write('MESH_DIMS = 2\n')
        f.write('MESH_TOP1 = TRUE\n')
        f.write('MESH_EPS1 = 1.0\n')
        f.write('MESH_RAD2 = 1\n')
        f.write('GEOM_FEAT = TRUE\n')
        f.write('VERBOSITY = 2')
    
    
    
    calc_dir = rpath+'/jigsaw/'
    
    #---------------------------------     
    logger.info('executing jigsaw\n')
    #--------------------------------- 
    
    #execute jigsaw
    ex=subprocess.Popen(args=['jigsaw {}'.format(tag+'.jig')], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
        
    with open(calc_dir+'err.log', 'w') as f: 
      for line in iter(ex.stderr.readline,b''): 
        f.write(line.decode(sys.stdout.encoding))   
#        logger.info(line.decode(sys.stdout.encoding))
    ex.stderr.close()            

    with open(calc_dir+'run.log', 'w') as f: 
      for line in iter(ex.stdout.readline,b''): 
        f.write(line.decode(sys.stdout.encoding))   
#        logger.info(line.decode(sys.stdout.encoding))
    ex.stdout.close()         
    
    #--------------------------------- 
    logger.info('Jigsaw FINISHED\n')
    #---------------------------------
   
    logger.info('..reading mesh\n')
    
    [nodes,edges,tria] = read_msh(rpath+'/jigsaw/'+ tag + '.msh')
    
    nodes = nodes.apply(pd.to_numeric)
    tria = tria.apply(pd.to_numeric)
    edges = edges.apply(pd.to_numeric)
    
    # Look for hanging nodes 
    tri3 = tria.values[:,:3]
    q = np.unique(tri3.flatten()) # all the unique nodes in elements

    dq = list(set(range(nodes.shape[0])) - set(q)) # the ones that are in gcrop but not in elems

    dq.sort()
    nodes = nodes.drop(dq) # drop nodes
    nodes = nodes.rename_axis('tag').reset_index() # reset index
    
    ### Re-index tessalation

    A, idxA = np.unique(nodes['tag'], return_inverse=True)
    B, idxB = np.unique(tria['a'], return_inverse=True)
    IDX = np.in1d(A,B)
    tria['a'] = idxA[IDX][idxB]
    B, idxB = np.unique(tria['b'], return_inverse=True)
    IDX = np.in1d(A,B)
    tria['b'] = idxA[IDX][idxB]
    B, idxB = np.unique(tria['c'], return_inverse=True)
    IDX = np.in1d(A,B)
    tria['c'] = idxA[IDX][idxB]
    
    # Drop invalid edges
    drop_e = edges.loc[edges.e1.isin(dq) | edges.e2.isin(dq)].index
    edges = edges.drop(drop_e).reset_index(drop=True)
    ### Re-index edges
    A, idxA = np.unique(nodes['tag'], return_inverse=True)
    B, idxB = np.unique(edges['e1'], return_inverse=True)
    IDX = np.in1d(A,B)
    edges['e1'] = idxA[IDX][idxB]
    B, idxB = np.unique(edges['e2'], return_inverse=True)
    IDX = np.in1d(A,B)
    edges['e2'] = idxA[IDX][idxB]
    #clean up
    nodes = nodes.drop('tag',axis=1)
   
       
    # Boundaries
    
    # LAND Boundaries (negative tag)    
    ib = 1
    isl = []
    for ik in range(edges.e3.min(),0):
        bb = np.unique(edges.loc[edges.e3 == ik,['e1','e2']].values.flatten())
        bf = pd.concat([nodes.loc[bb]],keys=['land_boundary_{}'.format(ib)])
        bf['idx'] = bb
        if ik >= bmindx: 
            bf['flag'] = 0
        else:
            bf['flag'] = 1

        if not bf.empty:
            isl.append(bf)
            ib += 1
    
    if isl :   
         
        landb = pd.concat(isl)
    
        land_i = sorted(landb.index.levels[0], key=lambda x: int(x.split('_')[2]))
    
    else:
        
        land_i = []
        
    # WATER Boundaries (positive tag)
    
    wb = 1
    wbs = []
    for ik in range(1,edges.e3.max()+1):
        bb = np.unique(edges.loc[edges.e3 == ik,['e1','e2']].values.flatten())
        bf = pd.concat([nodes.loc[bb]],keys=['open_boundary_{}'.format(wb)])
        bf['idx'] = bb

        wbs.append(bf)
        wb += 1
    
    if wbs:   
         
        openb = pd.concat(wbs)
        
        open_i = sorted(openb.index.levels[0], key=lambda x: int(x.split('_')[2]))    
    
        if openb.index.levels[0].shape[0] == 1: # sort the nodes if open box 
    
            pts = openb[['lon','lat']].values

            origin = [openb.mean()['lon'], openb.mean()['lat']]

            refvec = [1, 0]

            sps = sorted(pts, key=lambda po: clockwiseangle_and_distance(po, origin, refvec))

            sps = np.array(sps)
    
            # reorder openb
            fs = [tuple(lst) for lst in sps]
            ss = [tuple(lst) for lst in pts]
    
            index_dict = dict((value, idx) for idx,value in enumerate(ss))
            idx = [index_dict[x] for x in fs]
    
            openb = openb.iloc[idx]
    
    else:
        
        open_i = []
    
    #MAKE Dataset
    
    els = xr.DataArray(
          tria.loc[:,['a','b','c']].values,
          dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'], name='SCHISM_hgrid_face_nodes'
          )

    nod = nodes.loc[:,['lon','lat']].to_xarray().rename({'index':'nSCHISM_hgrid_node','lon':'SCHISM_hgrid_node_x','lat':'SCHISM_hgrid_node_y'})
    nod = nod.drop_vars('nSCHISM_hgrid_node')

    dep = xr.Dataset({'depth': (['nSCHISM_hgrid_node'], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    #open boundaries
    op=[]
    nop=[]
    o_label=[]
    for line in open_i:
        op.append(openb.loc[line,'idx'].values)
        nop.append(openb.loc[line,'idx'].size)
        o_label.append(line)

    xob = pd.DataFrame(op).T
    xob.columns = o_label

    oattr = pd.DataFrame({'label':o_label,'nps':nop})
    oattr['type'] = np.nan
    oattr.set_index('label', inplace=True, drop=True)

    #land boundaries

    lp=[]
    nlp=[]
    l_label=[]
    itype = []
    for line in land_i:
        lp.append(landb.loc[line,'idx'].values)
        nlp.append(landb.loc[line,'idx'].size)
        itype.append(landb.loc[line,'flag'].values[0])
        l_label.append(line)
    
    xlb = pd.DataFrame(lp).T
    xlb.columns = l_label

    lattr = pd.DataFrame({'label':l_label,'nps':nlp})
    lattr['type'] = itype
    lattr.set_index('label', inplace=True, drop=True)


    gr = xr.merge([nod,dep,els,xob.to_xarray(), xlb.to_xarray(), lattr.to_xarray(), oattr.to_xarray()]) # total
    
    
    logger.info('..done creating mesh\n')
    
    
    return gr
    
    