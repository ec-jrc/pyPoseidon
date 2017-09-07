import numpy as np
import pickle
import os
from Poseidon.model.grid import *
from Poseidon.model.dep import *
from Poseidon.utils.vis import *
import pyresample
import pandas as pd
import datetime
from netCDF4 import Dataset
import subprocess
import tempfile
import time

class data:   
       
                                    
    def animaker(self,var,vmax,vmin,title,label,units,step):
              
        flist=[]
       
        k=0

        for folder in self.folders:
                   
          # read array
          dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
          
          h = dat.variables[var][:step,1:-1,1:-1]
          ha = np.transpose(h,axes=(0,2,1)) #transpose lat/lon
          ww = np.broadcast_to(self.w == True, ha.shape) # expand the mask in time dimension
          z = np.ma.masked_where(ww==True,ha) #mask land
          
          a = anim(self.xh,self.yh,z,title=title,label=label,units=units,vrange=[vmin,vmax])
          
          with tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) as temp: # open temp file
              a.save(temp.name, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])      
              flist.append(temp.name)
          
          k+=1
          
        #save list 
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False) as temp:
              listname=temp.name
              for name in flist: temp.write('file    ' + name)
         
        #merge clips   
        outfile = tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) # open temp file
        ex=subprocess.Popen(args=['ffmpeg -y -f concat -i {} -c copy {}'.format(listname,outfile.name)], cwd=tempfile.gettempdir(), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
     #       out = video(temp.name, 'mp4')

        time.sleep(2)

        #cleanup
        for name in flist:
            os.remove(name)
        os.remove(listname)
        
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
        timef = []
        for folder in self.folders:
            with open(folder+'/info.pkl', 'r') as f:
                              info=pickle.load(f)
            
            dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
            vmax = max(dat.variables[var][:step].max(), vmax)
            vmin = min(dat.variables[var][:step].min(), vmin)
            t = dat.variables['time'][:step]
            for it in range(t.shape[0]):
                timef.append(info['date']+datetime.timedelta(seconds=np.int(t[it])))
            title = dat.variables[var].long_name
            units = dat.variables[var].units
            
        self.time = timef
         
        return self.animaker(var,vmax,vmin,title,timef,units,step)
        
    
   # def frame(self,slot):
        
                
            
    def point(self,plat,plon,**kwargs):
        
        var = kwargs.get('var', 'S1')
        
        plat=float(plat)
        plon=float(plon)
        
        j=np.abs(self.xh.data[0,:]-plon).argmin()
        i=np.abs(self.yh.data[:,0]-plat).argmin()
        
        #print i,j
        
        
        for folder in self.folders:
            
            with open(folder+'/info.pkl', 'r') as f:
                              info=pickle.load(f)
                   
            dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
        
            tv=dat.variables['time'][:]
            
            t = [info['date']+datetime.timedelta(seconds=np.int(it)) for it in tv]
            
            ndt = tv.shape[0]
        
            xb, yb = self.xh[i-3:i+4,j-3:j+4],self.yh[i-3:i+4,j-3:j+4] # retrieve nearby grid values
            orig = pyresample.geometry.SwathDefinition(lons=xb,lats=yb) # create original swath grid
            
            
            targ = pyresample.geometry.SwathDefinition(lons=np.array([plon,plon]),lats=np.array([plat,plat])) #  point
            
            pval=[]
            for it in range(ndt):
                # note the transposition of the D3D variables
                vals = np.ma.masked_array(dat.variables[var][it,j-2:j+5,i-2:i+5],self.w[i-3:i+4,j-3:j+4]) # values as nearby points ! land in masked !!
                #print vals
                vals.fill_value=999999
                
                s = pyresample.kd_tree.resample_nearest(orig,vals,targ,radius_of_influence=50000,fill_value=999999)
                pval.append(s[0])
                
            return pd.DataFrame({'time':t,var:pval})    
        
        