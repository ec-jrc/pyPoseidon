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
import folium
from pyPoseidon.utils.get_value import get_value
import xarray
import shutil
import holoviews as hv
import geoviews as gv
from cartopy import crs
import xarray as xr




FFWriter = animation.FFMpegWriter(fps=30, extra_args=['-vcodec', 'libx264','-pix_fmt','yuv420p'])
   
class data:   
       
                                    
    def animaker(self,var,vmax,vmin,title,units,step,savepath=None):
              
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
              a.save(temp.name, writer = FFWriter)      
              flist.append(temp.name)
          
          k+=1
          
        self.time = timef
        #save list 
        with tempfile.NamedTemporaryFile(mode='w+t', suffix='.txt', delete=False) as temp:
              listname=temp.name
              for name in flist: temp.write('file   {}\n'.format(name))
         
        #merge clips   
        outfile = tempfile.NamedTemporaryFile(suffix='.mp4', delete=False) # open temp file
        ex=subprocess.Popen(args=['ffmpeg -y -f concat -safe 0 -i {} -c copy {}'.format(listname,outfile.name)], cwd=tempfile.gettempdir(), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
     #       out = video(temp.name, 'mp4')

        time.sleep(5)

        #cleanup
        for name in flist:
            os.remove(name)
        os.remove(listname)

        plt.ion()
#       print outfile.name
#       print os.path.isfile(outfile.name) 
#       print os.path.getsize(outfile.name)
        if savepath:
           shutil.move(outfile.name,savepath)
           return
        else:
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
        
        self.dx = self.xh.data[0,1]-self.xh.data[0,0]
        self.dy = self.yh.data[1,0]-self.yh.data[0,0]        
         
        self.variables = dat.variables.keys()  
                 
        self.ds = xarray.open_mfdataset(nfiles)
                
    def hview(self,var,**kwargs):
        
        vart = np.squeeze(self.ds[var].values)
        vart = np.transpose(vart,axes=(0,2,1))

        lon = self.ds.XZ.values[1:-1,1:-1].T[0,:]

        lat = self.ds.YZ.values[1:-1,1:-1].T[:,0]

        times=self.ds.time.values

        ds=xr.Dataset({var:(['time','latitude','longitude'], vart[:,1:-1,1:-1])},
                     coords={'longitude':(lon),'latitude':(lat),'time':times})

        return hv.Dataset(ds,kdims=['time','longitude','latitude'],vdims=[var])
   
   
    def gview(self,var,**kwargs):
        
        vart = np.squeeze(self.ds[var].values)
        vart = np.transpose(vart,axes=(0,2,1))

        lon = self.ds.XZ.values[1:-1,1:-1].T[0,:]

        lat = self.ds.YZ.values[1:-1,1:-1].T[:,0]

        times=self.ds.time.values

        ds=xr.Dataset({var:(['time','latitude','longitude'], vart[:,1:-1,1:-1])},
                     coords={'longitude':(lon),'latitude':(lat),'time':times})

        return gv.Dataset(ds,kdims=['time','longitude','latitude'],vdims=[var])
    
                                   
    def movie(self,var,**kwargs):
         
        step = kwargs.get('step', None)
        savepath = kwargs.get('savepath', None)
         
        vmax=0.
        vmin=0. 
        for folder in self.folders:

            dat=Dataset(folder+'/'+'trim-'+self.info['tag']+'.nc')
            vmax = max(dat.variables[var][:step].max(), vmax)
            vmin = min(dat.variables[var][:step].min(), vmin)
            title = dat.variables[var].long_name
            units = dat.variables[var].units
            
         
        return self.animaker(var,vmax,vmin,title,units,step,savepath=savepath)
        
    
   # def frame(self,slot):
   
    def maxval(self,var,**kwargs):
       
       step = kwargs.get('step', None)
       
       vmax = np.amax(np.transpose(otp.ds.S1.values,axes=(0,2,1)), axis=0)           
   
       return np.ma.masked_array(vmax,self.xh.mask)
       
               
class point:
    
    def __init__(self,**kwargs):
                        
        self.lon = kwargs.get('lon', None) 
        self.lat = kwargs.get('lat', None) 
        self.data = kwargs.get('data', None)
        self.name = kwargs.get('name', None)
            
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
                elif method == 'gauss':
                    s = pyresample.kd_tree.resample_gauss(orig,vals,targ,radius_of_influence=50000,fill_value=999999,sigmas=25000)
                
                pval.append(s[0])
                
                
            tdat = pd.DataFrame({'time':t,var:pval})
                        
            frames.append(tdat)


        pdata = pd.concat(frames)
        setattr(self, var, pdata.set_index(['time']))
        
    def on_map(self,**kwargs):
        
            from folium.features import DivIcon
        
            name = get_value(self,kwargs,'name',None)
            marker = kwargs.get('marker', None)
                        
            if marker == None:
                attrs = folium.Marker([self.lat,self.lon],popup=name)
            elif (marker == 'circle'):
                attrs = folium.CircleMarker(location=[self.lat,self.lon],
                                radius = 5,
                                fill=True,
                                color='black',
                                popup=name) 
        
            return attrs
        
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
          
    def on_map(self,**kwargs):
        
        from folium.features import DivIcon
        
        dyy = self.data.dy/2 # staggered shift
        dxx = self.data.dx/2
        
        dfp=pd.DataFrame([{'ID':self.name,'kind':'pressure point','i':self.i,'j':self.j,'lon':self.ilon,'lat':self.ilat}])

        dfp = dfp.set_index('ID')

        htmlp = dfp.to_html(classes='table table-striped table-hover table-condensed table-responsive')#,index=None)
        
        popup_p = folium.Popup(htmlp)
                
        attrs=[]
        
        cell=folium.features.RectangleMarker(bounds=[[self.ilat-dyy, self.ilon-dxx], [self.ilat+dyy, self.ilon+dxx]], \
                                color='black', fill_color='black', fill_opacity=.4, popup='cell')
        attrs.append(cell)
        
        pp = folium.CircleMarker(location=[self.ilat,self.ilon],
                            radius = 6,
                            fill=True,
                            color='blue',
                            popup=popup_p)
                            
        attrs.append(pp)
        
 
        dfg=pd.DataFrame([{'ID':self.name, 'kind':'grid point','i':self.i,'j':self.j,'lon':self.ilon+dxx,'lat':self.ilat+dyy}])

        dfg = dfg.set_index('ID')
 
        htmlg = dfg.to_html(classes='table table-striped table-hover table-condensed table-responsive')#,index=self.name)
 
        popup_g = folium.Popup(htmlg)
 
                
        gp = folium.CircleMarker(location=[self.ilat+dyy,self.ilon+dxx],
                            radius = 3,
                            fill=True,
                            color='black',
                            popup=popup_g) 
                            
        attrs.append(gp)
        
        
        return attrs
        
                
#        folium.map.Marker(
#             [self.lat,self.lon],
#             icon=DivIcon(
#                 icon_size=(150,36),
#                 icon_anchor=(0,0),
#                 html='<div style="font-size: 12pt">Pressure point</div>',
#                 )
#             ).add_to(mapa)
         
#        folium.map.Marker(
#             [self.lat+self.data.dy,self.lon+self.data.dx],
#             icon=DivIcon(
#                 icon_size=(150,12),
#                 icon_anchor=(0,0),
#                 html='<div style="font-size: 12pt">Grid point</div>',
#                 )
#             ).add_to(mapa)
             
             
                            
        return mapa    
        
        
