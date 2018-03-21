import numpy as np
import datetime
import xarray as xr
import pandas as pd


class grid:
    
    def __init__(self, type=None, **kwargs):
        impl=None
        if type == 'r2d':
            self.impl = r2d(**kwargs)
        elif type == 'gsh':
            self.impl = gsh(**kwargs)            

class r2d(grid):
    """Regular 2d grid for d3d
    """
    def __init__(self, **kwargs):
        
        gx    = kwargs.get('x', None)
        gy    = kwargs.get('y', None)
        attrs = kwargs.get('attrs', {'Coordinate System': 'Spherical', 'alfori': 0.0, 'xori': 0.0, 'yori': 0.0})
        
        g = xr.Dataset({'lons': (['x', 'y'], gx),   
                        'lats': (['x', 'y'], gy)})
        g.attrs = attrs
        
        self.grid = g
        

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
          
        g = xr.Dataset({'lons': (['x', 'y'], lons),   
                        'lats': (['x', 'y'], lats)})
                        
        g.attrs = {'Coordinate System': cs, 'alfori': alfori, 'xori': xori, 'yori': yori}
        
        return g
            
                    
    def to_file(self, filename, **kwargs):
        with open(filename,'w') as f:
            f.write('Coordinate System= {}\n'.format(self.grid.attrs['Coordinate System']))
            f.write('{} {}\n'.format(self.grid.lons.shape[1],self.grid.lons.shape[0]))
            f.write('{} {} {}\n'.format(self.grid.attrs['xori'],self.grid.attrs['yori'],self.grid.attrs['alfori']))
            for i in range(self.grid.lons.shape[0]):
                f.write('ETA=  {} '.format(i+1))
                f.write(' '.join(map(str, self.grid.lons[i,:].values)))
                f.write('\n')
            for i in range(self.grid.lats.shape[0]):
                f.write('ETA=  {} '.format(i+1))
                f.write(' '.join(map(str, self.grid.lats[i,:].values)))
                f.write('\n')

        


class gsh(grid):
    """Unstructured triangular grid from Jigsaw
    """
