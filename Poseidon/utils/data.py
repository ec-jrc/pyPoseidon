import numpy as np
import os

class data:
    
    def __init__(self,loc,**kwargs):
            
        folders = [os.path.abspath(name) for name in os.listdir(loc) if os.path.isdir(name)]
               
        if len(folders)==0:
            self.anim = anim(self,loc,**kwargs)
        else:    
            self.anim = anim_(self,folders,**kwargs)
                  
    def anim(self,folder,**kwargs):
        
        var = kwargs.get('var', 'S1')
           
        with open(folder+'/info.pkl', 'r') as f:
                     info=pickle.load(f) 
         
        tag=info['tag']
                 
        data=DataFile(folder+'/trim-'+tag+'.nc')   
        grid=Grid.fromfile(folder+'/trim-'+tag+'.grd')
        deb=Dep.read(folder+'/trim-'+tag+'.dep',grid.shape)
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
               
                   
    def anim_(self,folders,**kwargs):
       
        list=[]
       
        k=0 
              
        for folder in folders:
           
          a = anim(self,folder)
          
          filename = '/tmp/anim{}.mp4'.format(k)
          
          a.save(filename, fps=10, extra_args=['-vcodec','libx264','-pix_fmt','yuv420p'])
          
          list.append(filename)
          
        #save list 
         
        #merge clips   
        ex=subprocess.Popen(args=['ffmpeg -f concat -i {} -c copy /tmp/movie.mp4'.format(flist)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
     
     
        return '/tmp/movie.mp4'