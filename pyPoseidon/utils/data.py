import numpy as np
import pickle
import os
from pyPoseidon.model.grid import *
from pyPoseidon.model.dep import *
from pyPoseidon.model import mdf
from pyPoseidon.utils.vis import *
import pyresample
from pyresample import bilinear
import pandas as pd
import datetime
from netCDF4 import Dataset
import subprocess
import tempfile
import time

   
class data:   
       
                                    
    def animaker(self,var,vmax,vmin,title,units,step):
              
        flist=[]
       
        k=0

        timef = []
        
        plt.ioff()

        for folder in self.folders:

          # read netcdf
          dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
                   
          #read reference date
          inp, order = mdf.read(folder+'/'+self.info['tag']+'.mdf')
          refdat = inp['Itdate']
          refdate = datetime.datetime.strptime(refdat,'%Y-%m-%d')
          t = dat.variables['time'][:step]
          label = []
          for it in range(t.shape[0]):
                label.append(refdate+datetime.timedelta(seconds=np.int(t[it])))
                timef.append(refdate+datetime.timedelta(seconds=np.int(t[it])))
          
          h = dat.variables[var][:step,1:-1,1:-1]
          ha = np.transpose(h,axes=(0,2,1)) #transpose lat/lon
          ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
          z = np.ma.masked_where(ww==True,ha) #mask land
          
          a = anim(self.xh,self.yh,z,title=title,label=label,units=units,vrange=[vmin,vmax])
          
          with tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) as temp: # open temp file
              a.save(temp.name, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])      
              flist.append(temp.name)
          
          k+=1
          
        self.time = timef
        #save list 
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False) as temp:
              listname=temp.name
              for name in flist: temp.write('file   {}\n'.format(name))
         
        #merge clips   
        outfile = tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) # open temp file
        ex=subprocess.Popen(args=['ffmpeg -y -f concat -i {} -c copy {}'.format(listname,outfile.name)], cwd=tempfile.gettempdir(), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
     #       out = video(temp.name, 'mp4')

        time.sleep(5)

        #cleanup
        for name in flist:
            os.remove(name)
        os.remove(listname)

        plt.ion()
  #      print os.path.isfile(outfile.name) 
  #      print os.path.getsize(outfile.name)
        
        return video(outfile.name, 'mp4')
        
    def __init__(self,folders,**kwargs):
        
        self.folders = folders #[os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]
               
        with open(self.folders[0]+'/info.pkl', 'r') as f:
                          self.info=pickle.load(f)  
                        
        grid=Grid.fromfile(self.folders[0]+'/'+self.info['tag']+'.grd')
            
        deb=Dep.read(self.folders[0]+'/'+self.info['tag']+'.dep',grid.shape)
        
        self.dem = deb
        
        self.grid = grid
            
        #create mask
            
        d=deb.val[1:-1,1:-1]
        self.w=np.isnan(d)
            
            # read array
        dat=Dataset(self.folders[0]+'/'+'trim-'+self.info['tag']+'.nc') 
            
        xz = dat.variables['XZ'][:]
        xz = xz.T[1:-1,1:-1]
        yz = dat.variables['YZ'][:]   
        yz = yz.T[1:-1,1:-1]  
            
            
            
        self.xh = np.ma.masked_array(xz, self.w) #mask land
        self.yh = np.ma.masked_array(yz, self.w)        
         
        self.variables = dat.variables.keys()   
                           
    def movie(self,var,**kwargs):
         
        step = kwargs.get('step', None)
         
        vmax=0.
        vmin=0. 
        for folder in self.folders:

            dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
            vmax = max(dat.variables[var][:step].max(), vmax)
            vmin = min(dat.variables[var][:step].min(), vmin)
            title = dat.variables[var].long_name
            units = dat.variables[var].units
            
         
        return self.animaker(var,vmax,vmin,title,units,step)
        
    
   # def frame(self,slot):
   
    def maxval(self,var,**kwargs):
       
       step = kwargs.get('step', None)
       
       dat=Dataset(self.folders[0]+'/'+'trim-'+self.info['tag']+'.nc')
       h = dat.variables[var][:step,1:-1,1:-1]
       ha = np.transpose(h,axes=(0,2,1)) #transpose lat/lon
       vmax = np.amax(ha,axis=0) # max values of all times

       
       for folder in self.folders[1:]:
           # read netcdf
           dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
           h = dat.variables[var][:step,1:-1,1:-1]
           ha = np.transpose(h,axes=(0,2,1)) #transpose lat/lon
           
           hmax = np.amax(ha,axis=0) # max values of all times
         
           vmax = np.maximum(hmax,vmax)
           
   
       return np.ma.masked_array(vmax,self.xh.mask)
       
       
    def get_data(self,var,**kwargs): 
          
       step = kwargs.get('step', None)
       
       
       stor = []
       for folder in self.folders:
           # read netcdf
           dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
           h = dat.variables[var][:step,1:-1,1:-1]
           ha = np.transpose(h,axes=(0,2,1)) #transpose lat/lon
           
           stor.append(ha)
           
       return  np.hstack(stor)  
           
        
class point:
    
    def __init__(self,**kwargs):
                        
        self.lon = kwargs.get('lon', None) 
        self.lat = kwargs.get('lat', None) 
        self.data = kwargs.get('data', None)
            
    def tseries(self,**kwargs):
        
        var = kwargs.get('var', None)
        step = kwargs.get('step', None)
        method = kwargs.get('method', 'nearest')
        
        plat=float(self.lat)
        plon=float(self.lon)
        
        j=np.abs(self.data.xh.data[0,:]-plon).argmin()
        i=np.abs(self.data.yh.data[:,0]-plat).argmin()
                
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
        
            xb, yb = self.data.xh[i-3:i+4,j-3:j+4],self.data.yh[i-3:i+4,j-3:j+4] # retrieve nearby grid values
            orig = pyresample.geometry.SwathDefinition(lons=xb,lats=yb) # create original swath grid
                                   
            targ = pyresample.geometry.SwathDefinition(lons=np.array([plon,plon]),lats=np.array([plat,plat])) #  point
            
            
            pval=[]
            for it in range(ndt):
                # note the transposition of the D3D variables
                svalues = dat.variables[var][it,j-2:j+5,i-2:i+5]
                mask = xb.mask
                vals = np.ma.masked_array(svalues,mask=mask) # values as nearby points ! land in masked !!
                #print vals
                vals.fill_value=999999
                
                if method == 'nearest':
                    s = pyresample.kd_tree.resample_nearest(orig,vals,targ,radius_of_influence=50000,fill_value=999999)
                elif methos == 'gauss':
                    s = pyresample.kd_tree.resample_gauss(orig,vals,targ,radius_of_influence=50000,fill_value=999999,sigmas=25000)
                
                pval.append(s[0])
                
                
            tdat = pd.DataFrame({'time':t,var:pval})
                        
            frames.append(tdat)


        pdata = pd.concat(frames)
        self.obs = pdata.set_index(['time'])
        
        
    def nearest_node(self,**kwargs):
                       
        plat=float(self.lat)
        plon=float(self.lon)
        
        j=np.abs(self.data.xh.data[0,:]-plon).argmin()
        i=np.abs(self.data.yh.data[:,0]-plat).argmin()
        
        self.i = i
        self.j = j
        
        self.ilon = self.data.xh.data[i,j]
        self.ilat = self.data.yh.data[i,j]
        
        
    def node_data(self,**kwargs):
        
        i = kwargs.get('i', self.i) 
        j = kwargs.get('j', self.j)  
        var = kwargs.get('var', None) 
        step = kwargs.get('step', None)   
        
        ilon, ilat = self.data.xh[i,j],self.data.yh[i,j] # retrieve nearby grid values
        
        p = point(lon=ilon,lat=ilat,data=self.data)
        
        
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
        p.obs = pdata.set_index(['time'])
        return p
    
    
    def next_node_data(self,**kwargs):
               
        var = kwargs.get('var', None)
        step = kwargs.get('step', None)
        nps = kwargs.get('npoints', 1)
        
        xx,yy =np.meshgrid(np.arange(-nps,nps+1),np.arange(-nps,nps+1))
        
        frames = []
        
        for l,k in zip((xx+self.i).flatten(),(yy+self.j).flatten()):
            p = self.node_data(i=l,j=k,var=var,step=step)
            frames.append(p)
     
        return frames 
          
        