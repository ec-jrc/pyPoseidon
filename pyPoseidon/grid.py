"""
Grid module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import datetime
import xarray as xr
import pandas as pd
import sys
from .jigsaw import *
from .ugmsh import *
import logging
import f90nml
import os
import subprocess
from pyPoseidon.utils.verify import *
import multiprocessing

NCORES = max(1, multiprocessing.cpu_count() - 1)

logger = logging.getLogger('pyPoseidon')


DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    


def grid(type=None, **kwargs):
    if type == 'r2d':
        return r2d(**kwargs)
    elif type == 'tri2d':
        return tri2d(**kwargs)            

class r2d():
    """Regular 2d grid for d3d
    """
    def __init__(self, **kwargs):
        
        grid_file  = kwargs.get('grid_file', None)
                    
        if grid_file :             
            
            self.Dataset = self.read_file(grid_file)
        
        else:
            
            geometry = kwargs.get('geometry', None)
    
            try:
                lon_min = geometry['lon_min']
                lon_max = geometry['lon_max']
                lat_min = geometry['lat_min']
                lat_max = geometry['lat_max']
            except:
                logger.error('geometry not set properly\n')
                sys.exit(1)
            
            resolution = kwargs.get('resolution', None)
            
            ni=int(round((lon_max-lon_min)/resolution)) #these are cell numbers
            nj=int(round((lat_max-lat_min)/resolution))
  
            lon_max=lon_min+ni*resolution #adjust max lon to much the grid
            lat_max=lat_min+nj*resolution

            # set the grid 
            x=np.linspace(lon_min,lon_max,ni)
            y=np.linspace(lat_min,lat_max,nj)
            gx,gy=np.meshgrid(x,y)

            attrs = kwargs.get('attrs', {'Coordinate System': 'Spherical', 'alfori': 0.0, 'xori': 0.0, 'yori': 0.0})
        
            g = xr.Dataset({'lons': (['y', 'x'], gx),   
                        'lats': (['y', 'x'], gy)},
                        coords={'x': ('x', gx[0,:]),   
                                'y': ('y', gy[:,0])})         
            
            g.attrs = attrs
        
            self.Dataset = g
        

    @staticmethod
    def read_file(filename, **kwargs):
        
        logger.info('read grid file {}'.format(filename))
        
        header=pd.read_csv(filename,nrows=3,header=None,comment='*')
        cs = header.loc[0,0].split('=')[1].strip()
        ni,nj = header.loc[1,0].split(' ')
        ni,nj = int(ni),int(nj)
        alfori, xori, yori = header.loc[2,0].split(' ')
        
        d = pd.read_csv(filename,header=2,comment='*',delim_whitespace=True,engine='python',na_values='ETA=')
        d = d.reset_index()
        data = d.values[~np.isnan(d.values)]
        data=np.array(data)
        data = data.reshape(2,nj,ni+1) # including the row index
        #clean up the row index
        data = data[:,:,1:]
        
        lons=data[0,:,:]
        lats=data[1,:,:]      
                        
        g = xr.Dataset({'lons': (['y', 'x'], lons),   
                        'lats': (['y', 'x'], lats)},
                         coords={'x': ('x', lons[0,:]),   
                                 'y': ('y', lats[:,0])})         
                        
                        
        g.attrs = {'Coordinate System': cs, 'alfori': alfori, 'xori': xori, 'yori': yori}
        
        return g
            
                    
    def to_file(self, filename, **kwargs):
        
        logger.info('writing grid to file {}'.format(filename))
        
        with open(filename,'w') as f:
            f.write('Coordinate System= {}\n'.format(self.Dataset.attrs['Coordinate System']))
            f.write('{} {}\n'.format(self.Dataset.lons.shape[1],self.Dataset.lons.shape[0]))
            f.write('{} {} {}\n'.format(self.Dataset.attrs['xori'],self.Dataset.attrs['yori'],self.Dataset.attrs['alfori']))
            for i in range(self.Dataset.lons.shape[0]):
                f.write('ETA=  {} '.format(i+1))
                f.write(' '.join(map(str, self.Dataset.lons[i,:].values)))
                f.write('\n')
            for i in range(self.Dataset.lats.shape[0]):
                f.write('ETA=  {} '.format(i+1))
                f.write(' '.join(map(str, self.Dataset.lats[i,:].values)))
                f.write('\n')

        
    
class tri2d():
    """Unstructured triangular 2d grid
    """
    def __init__(self, **kwargs):
                    
        grid_file  = kwargs.get('grid_file', None)
        grid_generator = kwargs.get('grid_generator', 'jigsaw')
                    
        if grid_file: 
              
            self.Dataset = self.read_file(grid_file)
            
        elif grid_generator == 'gmsh':
    
            g = gmsh_(**kwargs) # create grid with JIGSAW
    
            self.Dataset = g
            
        elif grid_generator == 'jigsaw':
            
            g = jigsaw(**kwargs) # create grid with JIGSAW
    
            self.Dataset = g
            
        else:
            
            self.Dataset = None
            
            
     
     
    @staticmethod
    def read_file(hgrid,**kwargs):
                
        logger.info('read grid file {}'.format(hgrid))
                
        #read file
        df = pd.read_csv(hgrid, header=0, names=['data'], index_col=None, low_memory=False)
        
        #extract number of elements, number of nodes
        ni,nj = df.iloc[0].str.split()[0]
        ni=int(ni)
        nj=int(nj) 
               
        #read lon,lat,depth for all nodes
        q = pd.DataFrame(df.loc[1:nj,'data'].str.split().values.tolist())
        q = q.drop(q.columns[0], axis=1)
        q = q.apply(pd.to_numeric)
      #  q.reset_index(inplace=True, drop=True)
        q.columns = ['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']
        q.index.name = 'nSCHISM_hgrid_node'
    
        #create xarray of grid
        grid = q.loc[:,['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y']].to_xarray()
        grid = grid.drop_vars('nSCHISM_hgrid_node')
        
        #create xarray of depth
        depth = q.loc[:,'depth'].to_xarray()
        depth = depth.drop_vars('nSCHISM_hgrid_node')
            
        #read connectivity
        e = pd.DataFrame(df.loc[nj+1:nj+ni,'data'].str.split().values.tolist())
        e = e.drop(e.columns[0], axis=1)
        e = e.apply(pd.to_numeric)
     #   e.reset_index(inplace=True, drop=True)
        e.columns = ['nv','a','b','c']
        e.loc[:,['a','b','c']] = e.loc[:,['a','b','c']] - 1 # convert to python (index starts from 0)
        
#        if e.nv.max() < 4:
#            e['d']=0
        
        #create xarray of tessellation
        els = xr.DataArray(
              e.loc[:,['a','b','c']].values,
              dims=['nSCHISM_hgrid_face', 'nMaxSCHISM_hgrid_face_nodes'], name='SCHISM_hgrid_face_nodes'
              )
                  
        #Open boundaries
        n0 = df[df.data.str.contains('open boundaries')].index
        n0 = n0.values[0]
        nob = df.loc[n0,'data'].split('=')[0].strip()
        nob = int(nob)
        nobn = df.loc[n0 + 1,'data'].split('=')[0].strip()
        nobn = int(nobn)
        
        onodes=[]
        ottr = []
        idx=n0 + 2
        for nl in range(nob):
            nn = df.loc[idx,'data'].split('=')[0].strip()
            nn = int(nn)
            label = df.loc[idx,'data'].split('=')[1]
            label = label[label.index('open'):]
            ottr.append([nn, label])
            nodes = df.loc[idx+1:idx+nn,'data']
            onodes.append(nodes.astype(int).values)
            idx = idx + nn + 1
        
        oinfo = pd.DataFrame(ottr)
        try:
            oinfo.columns = ['nps','label']
            oinfo.label = oinfo.label.str.rstrip()
            oinfo.label = oinfo.label.str.replace(' ', '_')
            oinfo.set_index('label', inplace=True, drop=True)
            oinfo = oinfo.apply(pd.to_numeric)
        except:
            pass


        obnodes = [j for i in onodes for j in i]
        dfo = pd.DataFrame({'node':obnodes,'type':-1,'id':0},index=np.arange(len(obnodes)))
        idx=1
        for l in range(len(onodes)):
            dfo.loc[dfo.node.isin(onodes[l]),'id'] = idx
            idx += 1
                
        #Land boundaries
        n1 = df[df.data.str.contains('land boundaries')].index
        n1 = n1.values[0]
        
        nlb = df.loc[n1,'data'].split('=')[0].strip()
        nlb = int(nlb)
        
        nlbn = df.loc[n1 + 1,'data'].split('=')[0].strip()
        nlbn = int(nlbn)
        

        lnodes=[]
        attr = []
        idx=n1 + 2
        for nl in range(nlb):
            nn, etype = df.loc[idx,'data'].split('=')[0].strip().split(' ')
            nn = int(nn)
            etype = int(etype)
            label = df.loc[idx,'data'].split('=')[1]
            label = label[label.index('land'):]
            attr.append([nn, etype, label])
            nodes = df.loc[idx+1:idx+nn,'data']
            lnodes.append(nodes.astype(int).values)
            idx = idx + nn + 1
            
        linfo = pd.DataFrame(attr)
        try:
            linfo.columns = ['nps','type','label']
            linfo.label = linfo.label.str.rstrip()
            linfo.label = linfo.label.str.replace(' ', '_') + '_' + linfo.type.astype(str)
            linfo.set_index('label', inplace=True, drop=True)
            linfo = linfo.apply(pd.to_numeric)
        except:
            pass
            
            
        lbnodes = [j for i in lnodes for j in i]
        dfl = pd.DataFrame({'node':lbnodes,'type':0,'id':0},index=np.arange(len(lbnodes)))
        idx=-1
        for l in range(len(lnodes)):
            dfl.loc[dfl.node.isin(lnodes[l]),'id'] = idx
            dfl.loc[dfl.node.isin(lnodes[l]),'type'] = linfo.iloc[l].type
            idx -= 1
        
        bbs = pd.concat([dfo,dfl])
        bbs.index.name = 'bnodes'
        
        bbs.node = bbs.node - 1 # start_index = 0  
            
        # merge to one xarray DataSet
        g = xr.merge([grid,depth,els,bbs.to_xarray()])
                    
        g.attrs = {}
        
        return g
    
    def to_file(self, filename, **kwargs):
        
        logger.info('writing grid to file {}'.format(filename))
        
        
        nn = self.Dataset.SCHISM_hgrid_node_x.size
        n3e = self.Dataset.nSCHISM_hgrid_face.size        
                
        with open(filename,'w') as f:
            f.write('\t uniform.gr3\n')
            f.write('\t {} {}\n'.format(n3e,nn))
        
        q = self.Dataset[['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']].to_dataframe()
        
        q.index = np.arange(1, len(q) + 1)
        
        q.to_csv(filename,index=True, sep='\t', header=None,mode='a', float_format='%.10f', columns=['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth'])   
        
        e = pd.DataFrame(self.Dataset.SCHISM_hgrid_face_nodes.dropna(dim='nMaxSCHISM_hgrid_face_nodes').values,columns=['a','b','c'])
        
        e['nv'] = e.apply(lambda row: row.dropna().size, axis=1)
        
        e.index = np.arange(1, len(e) + 1)
        
        e = e.dropna(axis=1).astype(int)
        
        e.loc[:,['a','b','c']] = e.loc[:,['a','b','c']] + 1 # convert to fortran (index starts from 1)
        
        e.to_csv(filename,index=True, sep='\t', header=None, mode='a', columns=['nv','a','b','c'])           
        
        
        bs = self.Dataset[['node','id','type']].to_dataframe()
        
        # open boundaries
        number_of_open_boundaries = bs.id.max().astype(int)
        number_of_open_boundaries_nodes=bs.loc[bs.id>0].shape[0]

        if number_of_open_boundaries > 0:
            with open(filename, 'a') as f:


                f.write('{} = Number of open boundaries\n'.format(number_of_open_boundaries))
                f.write('{} = Total number of open boundary nodes\n'.format(number_of_open_boundaries_nodes))


                for i in range(1, number_of_open_boundaries + 1):
                    dat = bs.loc[bs.id==i,'node'] + 1 #fortran
                    f.write('{} = Number of nodes for open boundary {}\n'.format(dat.size,i))
                    dat.to_csv(f,index=None,header=False)

        else:
    
            with open(filename, 'a') as f:

                f.write('{} = Number of open boundaries\n'.format(0))
                f.write('{} = Total number of open boundary nodes\n'.format(0))                                

        # land boundaries                      

        number_of_land_boundaries = bs.id.min().astype(int)
        number_of_land_boundaries_nodes=bs.loc[bs.id<0].shape[0]
        
        if number_of_land_boundaries < 0:
            with open(filename, 'a') as f:


                f.write('{} = Number of land boundaries\n'.format(-number_of_land_boundaries))
                f.write('{} = Total number of land boundary nodes\n'.format(number_of_land_boundaries_nodes))

                ik=1
                for i in range(-1 , number_of_land_boundaries - 1, -1):
                    dat_ = bs.loc[bs.id==i]
                    dat = dat_.node + 1 #fortran
                    dat_t = dat_.type.iloc[0].astype(int)
            
                    f.write('{} {} = Number of nodes for land boundary {}\n'.format(dat.size,dat_t,ik))
                    dat.to_csv(f,index=None,header=False)
                    ik += 1

        else:
    
            with open(filename, 'a') as f:

                f.write('{} = Number of land boundaries\n'.format(0))
                f.write('{} = Total number of land boundary nodes\n'.format(0))

    
        
    def validate(self,**kwargs):      
        
        #--------------------------------------------------------------------- 
        logger.info('start grid validation\n')
        #--------------------------------------------------------------------- 
        
        path = kwargs.get('rpath','./') 
        
        if not os.path.exists(path):
            os.makedirs(path)
               
        # save bctides.in
        bs = self.Dataset[['node','id','type']].to_dataframe()
        # open boundaries
        number_of_open_boundaries = bs.id.max().astype(int)
        number_of_open_boundaries_nodes=bs.loc[bs.id>0].shape[0]
        
        with open(path + 'bctides.in', 'w') as f:
            f.write('Header\n')
            f.write('{} {}\n'.format(0, 40.)) #  ntip tip_dp
            f.write('{}\n'.format(0)) #nbfr
            f.write('{}\n'.format(number_of_open_boundaries)) #number of open boundaries
            for i in range(1, number_of_open_boundaries + 1):
                nnodes = bs.loc[bs.id==i,'node'].shape[0]
                f.write('{} {} {} {} {}\n'.format(nnodes,2,0,0,0)) # number of nodes on the open boundary segment j (corresponding to hgrid.gr3), B.C. flags for elevation, velocity, temperature, and salinity
                f.write('{}\n'.format(0)) # ethconst !constant elevation value for this segment    
       
       #save vgrid.in
        with open(path + 'vgrid.in', 'w') as f:
            f.write('{}\n'.format(2)) #ivcor (1: LSC2; 2: SZ)
            f.write('{} {} {}\n'.format(2,1,1.e6)) #nvrt(=Nz); kz (# of Z-levels); hs (transition depth between S and Z)
            f.write('Z levels\n') # Z levels !Z-levels in the lower portion
            f.write('{} {}\n'.format(1,-1.e6)) #!level index, z-coordinates, z-coordinate of the last Z-level must match -hs 
            f.write('S levels\n') # S-levels below 
            f.write('{} {} {}\n'.format(40.,1.,1.e-4))  #constants used in S-transformation: h_c, theta_b, theta_f
            f.write('{} {}\n'.format(1,-1.)) #first S-level (sigma-coordinate must be -1)
            f.write('{} {}\n'.format(2,0.)) # levels index, sigma-coordinate, last sigma-coordinate must be 0

       # save params.nml
        config_file = DATA_PATH + 'param_val.nml'
        params = f90nml.read(config_file)
        params.write(path + 'param.nml',force=True)            
                         
        
        #save hgrid.gr3
        self.to_file(filename= path+'hgrid.gr3')
                             
        
        try:
            bin_path = os.environ['SCHISM']
        except:
            bin_path = kwargs.get('epath',None)
                                
        if bin_path is None:
            #------------------------------------------------------------------------------ 
            logger.warning('Schism executable path (epath) not given -> using default \n')
            #------------------------------------------------------------------------------
            bin_path = 'schism'
              
        ncores = NCORES
        calc_dir = path
                            
            
        with open(calc_dir + 'launchSchism.sh', 'w') as f:
            f.write('exec={}\n'.format(bin_path))
            f.write('mkdir outputs\n')
            f.write('mpirun -N {} $exec\n'.format(ncores))   
                 
          #make the script executable
        execf = calc_dir+'launchSchism.sh'
        mode = os.stat(execf).st_mode
        mode |= (mode & 0o444) >> 2    # copy R bits to X
        os.chmod(execf, mode)

            # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./launchSchism.sh'], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)#, bufsize=1)
            
        out, err = ex.communicate()[:]
        
        if 'successfully' in str(out):
            
            #--------------------------------------------------------------------- 
            logger.info('grid is validated for SCHISM\n')
            #--------------------------------------------------------------------- 
            return True
        else:
            logger.debug(str(out))
            #--------------------------------------------------------------------- 
            logger.info('grid fails.. exiting \n')
            #--------------------------------------------------------------------- 
            return False



    def verify(self,**kwargs):      
        
        shp = kwargs.get('coastlines',None) 
        
        if shp is not None:
            r=verify(self,shp)
            return r
        else:
            logger.warning('No coastlines provided')