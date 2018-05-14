import numpy as np
import pickle
import os
from pyPoseidon.model import *
from pyPoseidon.utils.vis import *
from pyPoseidon.utils.obs import obs
from pyPoseidon.grid import *
import pyresample
import pandas as pd
import datetime
from netCDF4 import Dataset
from pyPoseidon.utils.get_value import get_value
import holoviews as hv
import geoviews as gv
from cartopy import crs
import xarray as xr
import pyresample
import glob
import sys


FFWriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264','-pix_fmt','yuv420p'])
   
class data:   
       
        
    def __init__(self,folders,**kwargs):
        
        self.folders = folders #[os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]
        
        #check if many tags present
        ifiles = glob.glob(self.folders[0]+'/*_info.pkl')
        
        if len(ifiles) > 1:
            #--------------------------------------------------------------------- 
              sys.stdout.flush()
              sys.stdout.write('\n')
              sys.stdout.write('more than one configuration, specify tag argument \n')
              sys.stdout.flush()
            #--------------------------------------------------------------------- 
                    
        tag = kwargs.get('tag', None)
        
        if tag :
            ifile = self.folders[0]+'/'+tag+'_info.pkl'
        else:
            ifile = ifiles[0]
        
        #--------------------------------------------------------------------- 
        sys.stdout.flush()
        sys.stdout.write('\n')
        sys.stdout.write('reading data based on {} \n'.format(ifile))
        sys.stdout.flush()
        #--------------------------------------------------------------------- 
                       
        with open(ifile, 'r') as f:
                          self.info=pickle.load(f)  
                        
        grid=r2d.read_file(self.folders[0]+'/'+self.info['tag']+'.grd')
            
        deb=np.loadtxt(self.folders[0]+'/'+self.info['tag']+'.dep')
        
        #create mask
            
        d=deb[1:-1,1:-1]
        self.w=d==-999.
        
        b=deb[:-1,:-1]
        b[b==-999.]=np.nan
        
        self.dem = xr.Dataset({'bathymetry': (['latitude', 'longitude'], b)},
                    coords={'longitude': ('longitude', grid.lons[0,:]),   
                            'latitude': ('latitude', grid.lats[:,0])})
        
        self.grid = grid
            
            
            # read array
        
        nfiles = [folder +'/'+'trim-'+self.info['tag']+'.nc' for folder in folders]
        
        ds = xr.open_mfdataset(nfiles)
                
        #grid points
        xg = ds.XCOR.values.T
        yg = ds.YCOR.values.T   
        self.x = xg
        self.y = yg
        
        self.dx = self.x[0,1]-self.x[0,0]
        self.dy = self.y[1,0]-self.y[0,0]        
        
        #pressure points       
        xz = self.x - self.dx/2.
        yz = self.y - self.dy/2.
        
        xz = xz[1:-1,1:-1]  
        yz = yz[1:-1,1:-1]  
        
        self.xh = np.ma.masked_array(xz, self.w) #mask land
        self.yh = np.ma.masked_array(yz, self.w)
        
        #velocity points
        xv = xz + self.dx/2.
        yu = yz + self.dy/2.
                
        self.times=ds.time.values
        
        ovars = kwargs.get('vars', ['S1','U1','V1','WINDU','WINDV','PATM'])
        
        keys = [d for d in ds.variables.keys() if d in ovars]
        
        zkeys = [d for d in keys if d not in ['U1','V1']]
        
        dic = {}  
                
        for key in zkeys:
            h =  ds[key][:,1:-1,1:-1]
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            dic.update({key :(['time','YZ','XZ'], np.ma.masked_where(ww==True,ha))}) #mask land
        
        if 'U1' in keys:        
            # Get water velocities - U1
            h =  ds['U1'][:,0,1:-1,1:-1] 
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            wu = np.ma.masked_where(ww==True,ha) #mask land
            dic.update({'U1' :(['time','YU','XU'], wu)}) 
        
            #interpolate on XZ, YZ
            orig = pyresample.geometry.SwathDefinition(lons=xz,lats=yu) # original points
            targ = pyresample.geometry.SwathDefinition(lons=xz,lats=yz) # target grid
            uz = []
            for k in range(wu.shape[0]):
                uz.append(pyresample.kd_tree.resample_nearest(orig,wu[k,:,:],targ,radius_of_influence=50000))
         
            dic.update({'U1Z' :(['time','YZ','XZ'], uz)})
        
        if 'V1' in keys:        
            # Get water velocities - V1
            h =  ds['V1'][:,0,1:-1,1:-1] 
            ha = h.transpose(h.dims[0],h.dims[2],h.dims[1]) #transpose lat/lon
            ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
            wv = np.ma.masked_where(ww==True,ha) #mask land
            dic.update({'V1' :(['time','YV','XV'], wv)}) 
        
            #interpolate on XZ, YZ
            orig = pyresample.geometry.SwathDefinition(lons=xv,lats=yz) # original points
            targ = pyresample.geometry.SwathDefinition(lons=xz,lats=yz) # target grid
            vz = []
            for k in range(wv.shape[0]):
                vz.append(pyresample.kd_tree.resample_nearest(orig,wv[k,:,:],targ,radius_of_influence=50000))
      
            dic.update({'V1Z' : (['time','YZ','XZ'], vz)})
       
                 
        self.vars=xr.Dataset(dic,                                              
                     coords={'XZ':(xz[0,:]),'YZ':(yz[:,0]),
                             'XU':(xz[0,:]),'YU':(yu[:,0]),
                             'XV':(xv[0,:]),'YV':(yz[:,0]),                   
                     'time':self.times})
                     
                     
        #clean duplicates
        self.vars = self.vars.sel(time=~self.vars.indexes['time'].duplicated())
                     
        dic = self.info.copy()   # start with x's keys and values
        dic.update(kwargs)    # modifies z with y's keys and values & returns None
        
        if 'sa_date' not in dic.keys():
            dic.update({'sa_date':self.vars.time.values[0]})
            
        if 'se_date' not in dic.keys():
            dic.update({'se_date':self.vars.time.values[-1]})
               
        self.obs = obs(**dic)
                
    def hview(self,var,**kwargs):
        
        return hv.Dataset(self.vars[var])
   
   
    def gview(self,var,**kwargs):
        
        return gv.Dataset(self.vars[var])
    
                                   
    def frames(self,var,**kwargs):
        
        if len(var) == 1 :  
            return contour(self.xh,self.yh,self.vars[var[0]],self.times,**kwargs)
        elif len(var) == 2:
            return quiver(self.xh,self.yh,self.vars[var[0]],self.vars[var[1]],self.times,**kwargs)
            
            
                   
class point:
    
    def __init__(self,**kwargs):
                        
        self.lon = kwargs.get('lon', None) 
        self.lat = kwargs.get('lat', None) 
        self.data = kwargs.get('data', None)
            
    def tseries(self,**kwargs):
        
        var = kwargs.get('var', 'S1')
        method = kwargs.get('method', 'nearest')
        
        plat=float(self.lat)
        plon=float(self.lon)
 
        i=np.abs(self.data.xh[0,:].data-plon).argmin()
        j=np.abs(self.data.yh[:,0].data-plat).argmin()
                       
        xb, yb = self.data.xh[j-5:j+5,i-5:i+5],self.data.yh[j-5:j+5,i-5:i+5]
                       
        vals = self.data.vars[var][:,j-5:j+5,i-5:i+5].values
                        
        orig = pyresample.geometry.SwathDefinition(lons=xb,lats=yb) # create original swath grid
                                   
        targ = pyresample.geometry.SwathDefinition(lons=np.array([plon,plon]),lats=np.array([plat,plat])) #  point
        
        svals = []
                
        if method == 'nearest':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_nearest(orig,vals[k,:,:],targ,radius_of_influence=50000,fill_value=999)
                svals.append(s[0])
        elif method == 'gauss':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_gauss(orig,vals[k,:,:],targ,radius_of_influence=50000,fill_value=999,sigmas=25000)                
                svals.append(s[0])
                

        pdata = pd.DataFrame({'time':self.data.vars[var].time, self.data.vars[var].name : svals})
        setattr(self, self.data.vars[var].name, pdata.set_index(['time']))
        
        
class node:
    
    
    def __init__(self,**kwargs):

        self.i = kwargs.get('i', None) 
        self.j = kwargs.get('j', None)                         
        self.lon = kwargs.get('lon', None) 
        self.lat = kwargs.get('lat', None) 
        self.data = kwargs.get('data', None)
        self.name = kwargs.get('name', None)
        self.ilon = kwargs.get('ilon', None) 
        self.ilat = kwargs.get('ilat', None) 
    
    
    def index(self,**kwargs):
                       
        plat=float(self.lat)
        plon=float(self.lon)
        
        j=np.abs(self.data.xh.data[0,:]-plon).argmin()
        i=np.abs(self.data.yh.data[:,0]-plat).argmin()
        
        self.i = i
        self.j = j
        
        self.ilon = self.data.xh.data[i,j]
        self.ilat = self.data.yh.data[i,j]
        
        
    def tseries(self,**kwargs):
        
        try:
           i = kwargs.get('i', self.i) 
           j = kwargs.get('j', self.j)
           self.i = i
           self.j = j
           self.ilon, self.ilat = self.data.xh.data[i,j],self.data.yh.data[i,j] # retrieve nearby grid values
        except:
           self.index()
           i = self.i
           j = self.j

        if i == None :
           self.index()
           i = self.i
           j = self.j

        var = kwargs.get('var', get_value(self,kwargs,'var',None)) 
        step = kwargs.get('step', get_value(self,kwargs,'step',None))   

        frames=[] 
        
        for folder in self.data.folders:
            
            #read reference date
            inp, order = mdf.read(folder+'/'+self.data.info['tag']+'.mdf')
            refdat = inp['Itdate']
            refdate = datetime.datetime.strptime(refdat,'%Y-%m-%d')
                   
            dat=Dataset(folder+'/'+'trim-'+self.data.info['tag']+'.nc')
        
            tv=dat.variables['time'][:step]
            
            t = [refdate+datetime.timedelta(seconds=np.int(it)) for it in tv]
            
            ndt = tv.shape[0]
                      
            pval=dat.variables[var][:step,j,i] # retrieve variable data                
                
            tdat = pd.DataFrame({'time':t,var:pval})
                        
            frames.append(tdat)


        pdata = pd.concat(frames)
        setattr(self, var, pdata.set_index(['time']))
        
    
    def nearby_node_data(self,**kwargs):
               
        var = kwargs.get('var', None)
        step = kwargs.get('step', None)
        nps = kwargs.get('npoints', 1)
        
        xx,yy =np.meshgrid(np.arange(-nps,nps+1),np.arange(-nps,nps+1))
        
        frames = []
        
        m=0
        for l,k in zip((xx+self.i).flatten(),(yy+self.j).flatten()):
            p = node(data=self.data)
            p.tseries(i=l,j=k,var=var,step=step)
            p.name = 'point {}'.format(m)
            getattr(p, var).columns=[p.name]
            frames.append(p)
            m+=1
     
        return frames 
          
