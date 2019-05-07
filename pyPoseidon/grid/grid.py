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
        
        #extract number of elements, number of nodes
        ne,nn = pd.read_csv(hgrid,header=0,skiprows=1,nrows=0,delim_whitespace=True)
    
        ne=int(ne)
        nn=int(nn)        
        #read lon,lat,depth for all nodes
        q=pd.read_csv(hgrid,skiprows=2,header=None,delim_whitespace=True,engine='python',nrows=nn,names=['id','x','y','z'])        
        q=q.set_index(['id'])
    
        #read connectivity
        e = pd.read_csv(hgrid,skiprows=nn+2,header=None,delim_whitespace=True,engine='python',nrows=ne,names=['id','nv','a','b','c'])
        e=e.set_index(['id'])
    
        #Number of open boundaries
        nob = pd.read_csv(hgrid,header=None,skiprows=ne+nn+2,nrows=1,delimiter='=')[0]

        if nob.values > 0:
            #Total number of open boundary nodes
            ton = pd.read_csv(hgrid,header=None,skiprows=ne+nn+3,nrows=1,delimiter='=')[0]

            ns = ne + nn + 2 + 2
            nnob = pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[0].values[0]
            label = ' '.join(pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[1].str.split()[0][-3:])
            op = pd.read_csv(hgrid,header=None,skiprows=ns+1,nrows=nnob,names=[label])
            for i in range(1, nob):
                ns = ns + 1 + nnob
                nnob = pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[0].values[0]
                label = ' '.join( pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[1].str.split()[0][-3:] )
                op = pd.concat([op,pd.read_csv(hgrid,header=None,skiprows=ns+1,nrows=nnob,names=[label])], axis=1)

            op.index.name = 'id' # set index name to match q,e
            op.index = op.index + 1 #  ''
        
        else:
            ns = ne + nn + 3
            nnob = 0
            op = pd.DataFrame({})
        
        #Number of land boundaries
        nlb = pd.read_csv(hgrid,header=None,skiprows=ns + 1 + nnob,nrows=1,delimiter='=')[0]

        if nlb.values > 0:
            #Total number of land boundary nodes
            tln = pd.read_csv(hgrid,header=None,skiprows=ns + 1 + nnob + 1,nrows=1,delimiter='=')[0]

            ns = ns + 1 + nnob + 2
            nnlb, btype = pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[0].str.split()[0]
            nnlb = int(nnlb)
            label = ' '.join(pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[1].str.split()[0][-3:])
            bt=[btype]
            lp = pd.read_csv(hgrid,header=None,skiprows=ns+1,nrows=nnlb,names=[label])
            for i in range(1, nlb.values[0]):
                ns = ns + 1 + nnlb
                nnlb, btype  = pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[0].str.split()[0]
                nnlb = int(nnlb)

                bt.append(btype)
                label = ' '.join( pd.read_csv(hgrid,header=None,skiprows=ns,nrows=1,delimiter='=')[1].str.split()[0][-3:] )
                lp = pd.concat([lp,pd.read_csv(hgrid,header=None,skiprows=ns+1,nrows=nnlb,names=[label])], axis=1)

            lp.index.name = 'id' # set index name to match q,e
            lp.index = lp.index + 1 #  ''
        
            #DataFrame with the type of land boundary
            bbt = pd.DataFrame({'type':bt})
            bbt.index.name = 'id'
            bbt.index = bbt.index + 1
        
        else:
            lp=pd.DataFrame({})
            bbt = pd.DataFrame({})
    
        # merge to one xarray DataSet
        g = xr.merge([q.to_xarray(), e.to_xarray(), op.to_xarray(), lp.to_xarray(), bbt.to_xarray()])
                    
        g.attrs = {}
    
        return g
    
    def to_file(self, filename, **kwargs):
        
        nn = self.Dataset.x[np.isfinite(self.Dataset.x.values)].size
        n3e = self.Dataset.a.size        
                
        with open(filename,'w') as f:
            f.write('\t uniform.gr3\n')
            f.write('\t {} {}\n'.format(n3e,nn))
        
        q = self.Dataset.to_dataframe().loc[:,['x','y','z']].dropna()
        
        q.to_csv(filename,index=True, sep='\t', header=None,mode='a', float_format='%.10f', columns=['x','y','z'])   
        
        e = self.Dataset.to_dataframe().loc[:,['nv','a','b', 'c']] 
            
        e.to_csv(filename,index=True, sep='\t', header=None, mode='a', columns=['nv','a','b','c'])           
        
        # open boundaries
        keys = [k for k in self.Dataset.variables.keys() if 'open' in k]

        if keys :
        
            obound = self.Dataset.to_dataframe().loc[:,keys] # get the dataframe

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

        keys = [k for k in self.Dataset.variables.keys() if 'land' in k]

        if keys :

            lbound = self.Dataset.to_dataframe().loc[:,keys] # get the dataframe

            nlb = lbound.shape[1] # number of boundaries

            lps = (~lbound.isna()).sum() # number of nodes for each boundary
        
            with open(filename, 'a') as f:
                f.write('{} = Number of land boundaries\n'.format(nlb))
                f.write('{} = Total number of land boundary nodes\n'.format(lps.sum()))
                for i in range(nlb):
                    dat = lbound['land boundary {}'.format(i + 1)].dropna().astype(int)
                    f.write('{} {} = Number of nodes for land boundary {}\n'.format(dat.size,self.Dataset.type.values[i],i + 1))
                    dat.to_csv(f,index=None)

        

