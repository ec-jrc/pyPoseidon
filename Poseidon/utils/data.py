import numpy as np
import pickle
import os
from Poseidon.model.grid import *
from Poseidon.model.dep import *
from Poseidon.utils import *

class data:   
       
                  
    def anim(self,folder):
        
        var = kwargs.get('var', 'S1')
           
        with open(folder+'/info.pkl', 'r') as f:
                     info=pickle.load(f) 
         
        tag=info['tag']
                 
        data=DataFile(folder+'/trim-'+tag+'.nc')   
        grid=Grid.fromfile(folder+'/'+tag+'.grd')
        deb=Dep.read(folder+'/'+tag+'.dep',grid.shape)
        d=deb.val[:-1,:-1]
        w=np.isnan(d)
        x = grid.x.data
        y = grid.y.data
        h = data[var][:,:-1,:-1]
       
        # transpose prorerly
        ha=np.transpose(h,axes=(0,2,1))
        ww = np.broadcast_to(w == True, ha.shape)
       
        z = np.ma.masked_where(ww==True,ha)
        
        self.movie = anim(x,y,z[1:,:,:],title='Storm Surge',label='m',vrange=[z.min(),z.max()])
               
                   
    def anim_(self,folders):
       
        list=[]
       
        print folders
        
        k=0 
              
        for folder in folders:
         
          print folder
           
          a = anim(self,self,folder)
          
          filename = '/tmp/anim{}.mp4'.format(k)
          
          a.save(filename, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])
          
          list.append(filename)
          
        #save list 
         
        #merge clips   
        ex=subprocess.Popen(args=['ffmpeg -f concat -i {} -c copy /tmp/movie.mp4'.format(flist)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
     
     
        return '/tmp/movie.mp4'
        
        
    def __init__(self,loc,**kwargs):
        
            Mydict={'one':self.anim,'multiple':self.anim_}
            
            folders = [os.path.join(os.path.abspath(loc),name) for name in os.listdir(loc) if os.path.isdir(os.path.join(loc,name))]
            print folders
               
            if len(folders)==0:
                option = 'one'
            else:    
                option = 'multiple'
            
            print option    

            self.animation=Mydict[option](folders)
        
        