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
import matplotlib
from shapely import geometry, ops
import geopandas as gp
import xarray as xr
import logging
import os
import shapely
import subprocess
import sys

import pyPoseidon.dem as pdem
from .limgrad import *
        
        
#logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(levelname)-8s %(asctime)s:%(name)s:%(message)s')

file_handler = logging.FileHandler('jigsaw.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

sformatter = logging.Formatter('%(levelname)-8s %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(sformatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)
        


def jig(path='.',tag='jigsaw'):    
    
    fjig = path + '/' + tag+'.jig'
    
    with open(fjig,'w') as f:
        f.write('GEOM_FILE ={}\n'.format(tag+'-geo.msh'))
        f.write('MESH_FILE ={}\n'.format(tag+'.msh'))
        f.write('HFUN_FILE ={}\n'.format(tag+'-hfun.msh'))
        f.write('HFUN_SCAL = ABSOLUTE\n')
        f.write('HFUN_HMAX = Inf\n')
        f.write('HFUN_HMIN = 0.0\n')
        f.write('MESH_DIMS = 2\n')
        f.write('MESH_TOP1 = TRUE\n')
        f.write('MESH_EPS1 = 1.0\n')
        f.write('MESH_RAD2=1\n')
        f.write('VERBOSITY = 2')


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
        ie = df.loc[line].shape[0]+edges.shape[0]
        edges = edges.append(pd.DataFrame([[i,i + 1,0] for i in range(i0,ie)]))
        edges.iloc[-1][1]=i0

    # write header
    with open(fgeo,'a') as f:
        f.write('EDGE2={}\n'.format(edges.shape[0]))
    
    with open(fgeo, 'a') as f:
        edges.to_csv(f, index=False, header=0, sep=';')
        

def hfun(dem, path='.', tag='jigsaw', resolution_min=.5, resolution_max=10., dhdx=.15, imax=100, **kwargs):
    
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
    triang = matplotlib.tri.Triangulation(points[:,0], points[:,1])
    tri = triang.triangles
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
    
    # world polygons TODO user input
    world = gp.read_file(gp.datasets.get_path('naturalearth_lowres'))
    world = gp.GeoDataFrame.from_file('/Users/brey/DATA/COASTLINES/naturalearth/ne_10m_land/ne_10m_land.shp')
    world = world.explode()
    
    minlon = kwargs.get('minlon', None)
    maxlon = kwargs.get('maxlon', None)
    minlat = kwargs.get('minlat', None)
    maxlat = kwargs.get('maxlat', None)
    
    block = world.cx[minlon:maxlon,minlat:maxlat] #mask based on lat/lon window
    
    #create a polygon of the lat/lon window
    grp=geometry.Polygon([(minlon,minlat),(minlon,maxlat),(maxlon,maxlat),(maxlon,minlat)])

    #create a LineString of the grid
    grl=geometry.LineString([(minlon,minlat),(minlon,maxlat),(maxlon,maxlat),(maxlon,minlat),(minlon,minlat)])

    
    g = block.unary_union.symmetric_difference(grp) # get the dif from the world

    t = gp.GeoDataFrame({'geometry':g}) # make geoDataFrame    
    t['length']=t['geometry'][:].length #get length
    t = t.sort_values(by='length', ascending=0) #sort
    t = t.reset_index(drop=True)
    b = t.iloc[0].geometry #get the largest 
    
    #Here we extract all boundaries for the geo file
    dic={}
    for l in range(len(b.boundary)):
        lon=[]
        lat=[]
        for f in np.linspace(0.,b.boundary[l].length,1000): 
            cp = b.boundary[l].interpolate(f)
            lon.append(cp.x)
            lat.append(cp.y)

        dic.update({'line{}'.format(l):{'lon':lon,'lat':lat}})
    
    dict_of_df = {k: pd.DataFrame(v) for k,v in dic.items()} #concat
    df = pd.concat(dict_of_df, axis=0)
    df['z']=0
    df = df.drop_duplicates() # drop the repeat value on closed boundaries
    
    #open water boundaries
    water = b.boundary[0] - (b.boundary[0] - grl)
    land = b.boundary - grl # land boundaries!!

    tag = kwargs.get('tag', 'jigsaw')
    rpath = kwargs.get('rpath', '.')
    
    logger.info('Creating JIGSAW files\n')
    
    if not os.path.exists(rpath):
            os.makedirs(rpath)
    
    path = rpath+'/jigsaw/'
    if not os.path.exists(path):
        os.makedirs(path)
    
    
    jig(path=path,tag=tag)
    geo(df,path=path,tag=tag)
    
    dem_file = kwargs.get('dem_file', None)
    
    dem_dic={'minlon':minlon, # lat/lon window
         'maxlon':maxlon,
         'minlat':minlat,
         'maxlat':maxlat}
    if dem_file:
        dem_dic.update({'dem_file':dem_file})
         
    dem = pdem.dem(**dem_dic)
      
    
    hfun(dem, path=path,tag=tag)
    
    calc_dir = rpath+'/jigsaw/'
    
    #---------------------------------     
    logger.info('executing jigsaw\n')
    #--------------------------------- 
    
    #execute jigsaw
    ex=subprocess.Popen(args=[os.environ['CONDA_PREFIX']+'/envs/pyPoseidon/bin/jigsaw {}'.format(tag+'.jig')], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
        
    with open(calc_dir+'err.log', 'w') as f: 
      for line in iter(ex.stderr.readline,b''): 
        f.write(line.decode(sys.stdout.encoding))   
        logger.info(line.decode(sys.stdout.encoding))
    ex.stderr.close()            

    with open(calc_dir+'run.log', 'w') as f: 
      for line in iter(ex.stdout.readline,b''): 
        f.write(line.decode(sys.stdout.encoding))   
        logger.info(line.decode(sys.stdout.encoding))
    ex.stdout.close()         
    
    #--------------------------------- 
    logger.info('Jigsaw FINISHED\n')
    #---------------------------------
   
    logger.info('..reading mesh\n')
    
    [nodes,edges,tria] = read_msh(rpath+'/jigsaw/'+ tag + '.msh')
    
    nodes = nodes.apply(pd.to_numeric)
    tria = tria.apply(pd.to_numeric)
    
#    dem_dic.update({'grid_x':nodes.lon.values, 'grid_y':nodes.lat.values})
    
#    dem = pdem.dem(**dem_dic)
   
    # All elements
    polygons=[]
    for i in range(tria.shape[0]):
        temp = nodes.loc[tria.loc[i,['a','b','c']],['lon','lat']]
        temp['Coordinates'] = list(zip(temp.lon, temp.lat))
        poly = shapely.geometry.Polygon(list(temp.Coordinates.values))
        polygons.append(poly)
        
    gpt = gp.GeoDataFrame(geometry=polygons) # all polygons
    
    # get the boundary of the union
    cu = gpt.cascaded_union
    bdic={}
    for i in range(len(cu.boundary)):
        xb=cu.boundary[i].xy[0]
        yb=cu.boundary[i].xy[1]        
        bg = pd.DataFrame({'lon':xb, 'lat':yb})
        bdic.update({i:bg})
        
    bnodes  = pd.concat(bdic, axis=0)
    
    #sort out open/land
    
    gdf = gp.GeoDataFrame(
        bnodes, geometry=gp.points_from_xy(bnodes.lon, bnodes.lat))

    
    nodes['pos']= nodes['lon'].map(str) +',' + nodes['lat'].map(str) #create (lon,lat) for sorting the index
    
    #open boundaries
    openb={}
    for i in range(len(water)):
        openb.update({'open_boundary_{}'.format(i+1):gdf.loc[gdf.intersects(water[i].buffer(.1))]})
        
    openb = pd.concat(openb,axis=0)
    openb.index = openb.index.droplevel(1)
    openb['pos']= openb['lon'].map(str) +',' + openb['lat'].map(str)
    openb = openb.drop('geometry', axis=1).drop_duplicates()
    
    #get index 
    wid=[]
    for l in range(openb.index.levels[0].shape[0]):
        hd=openb.loc['open_boundary_{}'.format(l+1)]
        jj =hd.merge(nodes, how='inner', on=['pos'],left_index=True).index.values
        wid.append(jj)    
    
    for l in range(len(wid)):
        openb.loc['open_boundary_{}'.format(l+1),'idx']=wid[l]
    openb['idx']=openb['idx'].astype(int)
    
    #land boundaries
    landb={}
    for i in range(len(land)):
        landb.update({'land_boundary_{}'.format(i+1):gdf.loc[gdf.intersects(land[i].buffer(.1))]})
    
    landb = pd.concat(landb,axis=0)
    landb.index = landb.index.droplevel(1)
    landb['pos']= landb['lon'].map(str) +',' + landb['lat'].map(str)
    landb = landb.drop('geometry', axis=1).drop_duplicates()
    
    #get index
    lid=[]
    for l in range(landb.index.levels[0].shape[0]):
        idx=[]
        hd=landb.loc['land_boundary_{}'.format(l+1)]
        jj =hd.merge(nodes, how='inner', on=['pos'],left_index=True).index.values
        lid.append(jj)    
        
    for l in range(len(lid)):
        landb.loc['land_boundary_{}'.format(l+1),'idx']=lid[l]
    landb['idx']=landb['idx'].astype(int)
    
    
    #check dual use
    v = openb.merge(landb, how='inner', on=['pos'],right_index=True).index.values  
    for line,value in v:
        openb = openb.drop(value, level=1)
    
    #make Dataset
    nod = nodes.loc[:,['lon','lat']].to_xarray().rename({'lon':'SCHISM_hgrid_node_x','lat':'SCHISM_hgrid_node_y'}) # nodes

    els = xr.DataArray(
          tria.loc[:,['a','b','c','d']].values + 1 , # for index starting at one
          dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'], name='SCHISM_hgrid_face_nodes'
          )

    nod = nodes.loc[:,['lon','lat']].to_xarray().rename({'index':'nSCHISM_hgrid_node','lon':'SCHISM_hgrid_node_x','lat':'SCHISM_hgrid_node_y'})
    nod = nod.drop('nSCHISM_hgrid_node')

    dep = xr.Dataset({'depth': (['nSCHISM_hgrid_node'], np.zeros(nodes.shape[0]))})

    #open boundaries
    op=[]
    nop=[]
    o_label=[]
    for i in range(len(openb.index.levels[0])):
        line = 'open_boundary_{}'.format(i+1)
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
    for i in range(len(landb.index.levels[0])):
        line = 'land_boundary_{}'.format(i+1)
        lp.append(landb.loc[line,'idx'].values)
        nlp.append(landb.loc[line,'idx'].size)
        l_label.append(line)
        
    xlb = pd.DataFrame(lp).T
    xlb.columns = l_label

    lattr = pd.DataFrame({'label':l_label,'nps':nlp})
    lattr['type'] = 1.
    lattr.set_index('label', inplace=True, drop=True)

    

    gr = xr.merge([nod,dep,els,xob, xlb, lattr, oattr]) # total
    
    return gr
    
    