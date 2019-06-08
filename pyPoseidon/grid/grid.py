"""
Grid module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon, a software written by George Breyiannis (JRC E.1)
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import datetime
import xarray as xr
import pandas as pd
import sys

class grid:
    
    def __init__(self, type=None, **kwargs):
        impl=None
        if type == 'r2d':
            self.impl = r2d(**kwargs)
        elif type == 'tri2d':
            self.impl = tri2d(**kwargs)            

class r2d(grid):
    """Regular 2d grid for d3d
    """
    def __init__(self, **kwargs):
        
        grid_file  = kwargs.get('grid_file', None)
                    
        if grid_file :             
            
            self.Dataset = self.read_file(grid_file)
        
        else:
            
            minlon = kwargs.get('minlon', None)
            maxlon = kwargs.get('maxlon', None)
            minlat = kwargs.get('minlat', None)
            maxlat = kwargs.get('maxlat', None)
            resolution = kwargs.get('resolution', None)
            
            ni=int(round((maxlon-minlon)/resolution)) #these are cell numbers
            nj=int(round((maxlat-minlat)/resolution))
  
            maxlon=minlon+ni*resolution #adjust max lon to much the grid
            maxlat=minlat+nj*resolution

            # set the grid 
            x=np.linspace(minlon,maxlon,ni)
            y=np.linspace(minlat,maxlat,nj)
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

        
    
class tri2d(grid):
    """Unstructured triangular 2d grid
    """
    def __init__(self, **kwargs):
                    
        grid_file  = kwargs.get('grid_file', None)
                    
        if grid_file: 
              
            self.Dataset = self.read_file(grid_file)
        
        else:
    
            g = None # create grid with JIGSAW
    
            self.Dataset = g
    
    
    @staticmethod
    def read_file(hgrid,**kwargs):
                
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('reading grid from {}\n'.format(hgrid))
        sys.stdout.flush()
        
        #read file
        df = pd.read_csv(hgrid, header=0, names=['data'], index_col=None, low_memory=False)
        
        #extract number of elements, number of nodes
        ni,nj = df.iloc[0].str.split()[0]
        ni=int(ni)
        nj=int(nj) 
               
        #read lon,lat,depth for all nodes
        q = df.loc[1:nj,'data'].str.split('\t', 3, expand=True)
        q = q.drop(q.columns[0], axis=1)
        q = q.apply(pd.to_numeric)
      #  q.reset_index(inplace=True, drop=True)
        q.columns = ['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']
        q.index.name = 'nSCHISM_hgrid_node'
    
        #create xarray of grid
        grid = q.loc[:,['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y']].to_xarray()
        grid = grid.drop('nSCHISM_hgrid_node')
        
        #create xarray of depth
        depth = q.loc[:,'depth'].to_xarray()
        depth = depth.drop('nSCHISM_hgrid_node')
            
        #read connectivity
        e = df.loc[nj+1:nj+ni,'data'].str.split('\t', 5, expand=True)
        e = e.drop(e.columns[0], axis=1)
        e = e.apply(pd.to_numeric)
     #   e.reset_index(inplace=True, drop=True)
        e.columns = ['nv','a','b','c']
        if e.nv.max() < 4:
            e['d']=np.nan
        
        #create xarray of tessellation
        els = xr.DataArray(
              e.loc[:,['a','b','c','d']].values,
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
            oinfo.label = oinfo.label.str.replace(' ', '_')
            oinfo.set_index('label', inplace=True, drop=True)
            oinfo = oinfo.apply(pd.to_numeric)
        except:
            pass

        ops = pd.DataFrame(onodes).T
        try:
            ops.columns = oinfo.index
        except:
            pass
        
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
            linfo.label = linfo.label.str.replace(' ', '_')
            linfo.set_index('label', inplace=True, drop=True)
            linfo = linfo.apply(pd.to_numeric)
        except:
            pass
            
        lps = pd.DataFrame(lnodes).T
        lps.columns = linfo.index
            
            
        # merge to one xarray DataSet
        g = one = xr.merge([grid,depth,els,ops,lps,oinfo,linfo])
                    
        g.attrs = {}
    
        return g
    
    def to_file(self, filename, **kwargs):
        
        nn = self.Dataset.SCHISM_hgrid_node_x.size
        n3e = self.Dataset.nSCHISM_hgrid_face.size        
                
        with open(filename,'w') as f:
            f.write('\t uniform.gr3\n')
            f.write('\t {} {}\n'.format(n3e,nn))
        
        q = self.Dataset[['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']].to_dataframe()
        
        q.to_csv(filename,index=True, sep='\t', header=None,mode='a', float_format='%.10f', columns=['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth'])   
        
        e = pd.DataFrame(self.Dataset.SCHISM_hgrid_face_nodes.values,columns=['a','b','c','d'])
        
        e['nv'] = e.apply(lambda row: row.dropna().size, axis=1)
            
        e.to_csv(filename,index=True, sep='\t', header=None, mode='a', columns=['nv','a','b','c'])           
        
        # open boundaries
        keys = [k for k in self.Dataset.keys() if 'open' in k]

        if keys :
        
            obound = self.Dataset[keys].to_dataframe() # get the dataframe

            nob = obound.shape[1] # number of boundaries

            ops = (~obound.isna()).sum() # number of nodes for each boundary
        
            with open(filename, 'a') as f:
                f.write('{} = Number of open boundaries\n'.format(nob))
                f.write('{} = Total number of open boundary nodes\n'.format(ops.sum()))
                for i in range(nob):
                    dat = obound['open boundary {}'.format(i + 1)].dropna().astype(int)
                    f.write('{} = Number of nodes for open boundary {}\n'.format(dat.size,i+1))
                    dat.to_csv(f,index=None)
                                

        # land boundaries                      

        keys = [k for k in self.Dataset.keys() if 'land' in k]

        if keys :

            lbound = self.Dataset[keys].to_dataframe() # get the dataframe

            nlb = lbound.shape[1] # number of boundaries

            lps = (~lbound.isna()).sum() # number of nodes for each boundary
        
            with open(filename, 'a') as f:
                f.write('{} = Number of land boundaries\n'.format(nlb))
                f.write('{} = Total number of land boundary nodes\n'.format(lps.sum()))
                for i in range(nlb):
                    dat = lbound['land boundary {}'.format(i + 1)].dropna().astype(int)
                    f.write('{} {} = Number of nodes for land boundary {}\n'.format(dat.size,self.Dataset.type.values[i],i + 1))
                    dat.to_csv(f,index=None)

        

