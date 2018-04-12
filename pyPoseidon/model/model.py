import os
import datetime
import numpy as np
import pickle
import xml.dom.minidom as md
from shutil import copy2
import subprocess
import sys
import pkg_resources
from bunch import Bunch
import json
from collections import OrderedDict
import pandas as pd

#local modules
from pyPoseidon.grid import *
from bnd import *
import pyPoseidon
from pyPoseidon.meteo import meteo
from pyPoseidon.dem import dem
from pyPoseidon.utils.get_value import get_value

#retrieve the module path
#DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc')
DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    

# strings to be used 
le=['A','B']

nm = ['Z', 'A']

class model:
    impl = None
    def __init__(self, solver=None, **kwargs):
        if solver == 'd3d':
            self.impl = d3d(**kwargs)
        elif solver == 'schism':
            self.impl = schism(**kwargs)            
        
    def set(self,**kwargs):
        self.impl.set(**kwargs)
        
    def uvp(self,**kwargs):
        self.impl.uvp(**kwargs)
                  
    def pickle(self,**kwargs):
        self.impl.pickle(**kwargs)
            
    def run(self,**kwargs):
        self.impl.run(**kwargs)
        
    def save(self,**kwargs):
        self.impl.save(**kwargs)
    
    def read(self,**kwargs):
        self.impl.read(**kwargs)

    def force(self,**kwargs):
        self.impl.force(**kwargs)

    def bath(self,**kwargs):
        self.impl.bath(**kwargs)
        
    def output(self,**kwargs):
        self.impl.output(**kwargs)
        
    def bc(self,**kwargs):
        self.impl.bc(**kwargs)
                                   
    def tidebc(self,**kwargs):
        self.impl.tidebc(**kwargs)
        
class d3d(model):
    
    def __init__(self,**kwargs):
                
        self.minlon = kwargs.get('minlon', None)
        self.maxlon = kwargs.get('maxlon', None)
        self.minlat = kwargs.get('minlat', None)
        self.maxlat = kwargs.get('maxlat', None)
               
        start_date = kwargs.get('start_date', None)
        self.start_date = pd.to_datetime(start_date)
        
        if 'time_frame' in kwargs:
            time_frame = kwargs.get('time_frame', None)
            self.end_date = self.start_date + pd.to_timedelta(time_frame)
        else:
            end_date = kwargs.get('end_date', None)
            self.end_date = pd.to_datetime(end_date)
        
        if not hasattr(self, 'date'): self.date = self.start_date
        
        if self.end_date == None : sys.exit(1)
        
        self.tag = kwargs.get('tag', 'd3d')
        self.resolution = kwargs.get('resolution', None)
        self.ft1 = kwargs.get('ft1', None)
        self.ft2 = kwargs.get('ft2', None)
        self.dft = kwargs.get('dft', 1)       
        self.tide = kwargs.get('tide', False)
        self.atm = kwargs.get('atm', True)

        for attr, value in kwargs.iteritems():
                if not hasattr(self, attr): setattr(self, attr, value)
            
        self.model = self.__dict__.copy()    
        self.model['solver'] = self.__class__.__name__    
                      
    def set(self,**kwargs):

        gx = get_value(self,kwargs,'x',None)
        gy = get_value(self,kwargs,'y',None)    
        mdf_file = kwargs.get('mdf', None)  
        Tstart = self.start_date.hour*60     
        Tstop = int((self.end_date - self.start_date).total_seconds()/60)

        step = get_value(self,kwargs,'step',0)
        rstep = get_value(self,kwargs,'rstep',0)
        dt = get_value(self,kwargs,'Dt',1)
        
                                       
        resmin=self.resolution*60
              
        # computei ni,nj / correct lat/lon

        if gx :
            
          self.x=gx
          self.y=gy
              
          nj,ni=self.x.shape
          self.minlon=self.x.min()
          self.maxlon=self.x.max()
          self.minlat=self.y.min()
          self.maxlat=self.y.max()
          
        else:

          ni=int(round((self.maxlon-self.minlon)/self.resolution)) #these are cell numbers
          nj=int(round((self.maxlat-self.minlat)/self.resolution))
  
          self.maxlon=self.minlon+ni*self.resolution #adjust max lon to much the grid
          self.maxlat=self.minlat+nj*self.resolution

        # set the grid 
          x=np.linspace(self.minlon,self.maxlon,ni)
          y=np.linspace(self.minlat,self.maxlat,nj)
          gx,gy=np.meshgrid(x,y)

       
    #      ni=ni+1 # transfrom to grid points
    #      nj=nj+1
        
        self.ni=ni
        self.nj=nj
 
        # Grid 
          
        self.grid=grid(type='r2d',x=gx, y=gy)
        
        
        #mdf
        if mdf_file :
            self.mdf = pd.read_csv(mdf_file,sep='=')
        else:
            self.mdf = pd.read_csv(DATA_PATH+'default.mdf',sep='=')
        
        self.mdf = self.mdf.set_index(self.mdf.columns[0]) # set index
        
        
        mdfidx = self.mdf.index.str.strip() # store the stripped names
            
        #define grid file
        self.mdf.loc[self.mdf.index.str.contains('Filcco')]='#{}#'.format(self.tag+'.grd')
  
        #define enc file
        self.mdf.loc[self.mdf.index.str.contains('Filgrd')]='#{}#'.format(self.tag+'.enc')
  
        #define dep file
        self.mdf.loc[self.mdf.index.str.contains('Fildep')]='#{}#'.format(self.tag+'.dep')
  
        #define obs file
        self.mdf.loc[self.mdf.index.str.contains('Filsta')]='##' #b.tag+'.obs'
   
        # adjust ni,nj
        self.mdf.loc[self.mdf.index.str.contains('MNKmax')]='{} {} {}'.format(self.ni+1,self.nj+1,1)  # add one like ddb
  
        # adjust iteration date
        self.mdf.loc[self.mdf.index.str.contains('Itdate')]='#{}#'.format(self.date.strftime(format='%Y-%m-%d'))
  
        #set time unit
        self.mdf.loc[self.mdf.index.str.contains('Tunit')]='#M#'

        #adjust iteration start
        self.mdf.loc[self.mdf.index.str.contains('Tstart')]=Tstart
  
        #adjust iteration stop
        self.mdf.loc[self.mdf.index.str.contains('Tstop')]=Tstop
  
        #adjust time step
        self.mdf.loc[self.mdf.index.str.contains('Dt')]=[1.]
  
        #adjust time for output
        self.mdf.loc[self.mdf.index.str.contains('Flmap')]='{} {} {}'.format(Tstart,step,Tstop)
        self.mdf.loc[self.mdf.index.str.contains('Flhis')]='{} {} {}'.format(Tstart,dt,Tstop)
        self.mdf.loc[self.mdf.index.str.contains('Flpp')]='0 0 0'
        self.mdf.loc[self.mdf.index.str.contains('Flrst')]=rstep
  
        #time interval to smooth the hydrodynamic boundary conditions
        self.mdf.loc[self.mdf.index.str.contains('Tlfsmo')]=0.

        if not self.atm: self.mdf.loc['Sub1'] = ' '

        # set tide only run
        if self.tide :
            self.mdf.loc[self.mdf.index.str.contains('Filbnd')]='#{}#'.format(self.tag+'.bnd')
            self.mdf.loc[self.mdf.index.str.contains('Filana')]='#{}#'.format(self.tag+'.bca')
 #           if 'Tidfor' not in order: order.append('Tidfor')
 #           inp['Tidfor']=[['M2','S2','N2','K2'], \
 #                       ['K1','O1','P1','Q1'], \
 #                         ['-----------']]
 
          # specify ini file
        # if 'Filic' not in order: order.append('Filic')
        # inp['Filic']=basename+'.ini'

          # netCDF output
        if not 'FlNcdf' in mdfidx: 
            self.mdf.reindex(self.mdf.index.values.tolist()+['FlNcdf '])
            
        self.mdf.loc['FlNcdf '] = '#map his#'

        # Check for any other mdf variable in input
        for key,val in kwargs.iteritems():
            if key in mdfidx: self.mdf.loc[self.mdf.index.str.contains(key)] = val
        

        # get bathymetry
        self.bath()

        # get boundaries
        self.bc()
                
        #get meteo
        if self.atm :  self.force()
        
        #get tide
        if self.tide : self.tidebc()
        

        #meteo
    def force(self,**kwargs):
        z = self.__dict__.copy()
        z.update(kwargs)

        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:     
            check=[os.path.exists(z['rpath']+f) for f in ['u.amu','v.amv','p.amp']]   
            if (np.any(check)==False) :              
                self.meteo = meteo(**z)
            elif 'meteo' in flag : 
                self.meteo = meteo(**z)        
            else:
                sys.stdout.write('skipping meteo files ..\n')
        else:
            self.meteo = meteo(**z)
        
        #dem
    def bath(self,**kwargs):
        z = self.__dict__.copy()        
        
        z['grid_x'] = self.grid.impl.grid.lons.values
        z['grid_y'] = self.grid.impl.grid.lats.values
        
        z.update(kwargs) 
        
        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:
            check=[os.path.exists(z['rpath']+f) for f in ['{}.dep'.format(z['tag'])]]   
            if (np.any(check)==False) :
                self.dem = dem(**z)  
            elif 'dem' in flag :
                self.dem = dem(**z)
            else:
                sys.stdout.write('skipping dem files ..\n')
        else:
            self.dem = dem(**z)
            
        
    def bc(self,**kwargs):
        #define boundaries
        z = self.__dict__.copy()        
        
        z['lons'] = self.grid.impl.grid.lons[0,:]
        z['lats'] = self.grid.impl.grid.lats[:,0]
        
        try:
            ba = -self.dem.impl.ival.astype(np.float)
      # ba[ba<0]=np.nan
            z['dem']=ba
            z['cn']=10
        
            z.update(kwargs) 
                
            self.bound = box(**z)
        
        except:
            sys.stdout.write('boundary files not set..\n')
            

    def tidebc(self,**kwargs):
    
        self.tide = tide()
        for key,val in self.bound.__dict__.iteritems():
        
        # compute tide constituents
            tval = []
            if len(val) > 0. :                   
                blons=[]
                blats=[]
                for l1,l2 in val:
                    blons.append(self.grid.impl.grid.lons[l1[1]-1,l1[0]-1])   
                    blats.append(self.grid.impl.grid.lats[l1[1]-1,l1[0]-1])
                    blons.append(self.grid.impl.grid.lons[l2[1]-1,l2[0]-1])   
                    blats.append(self.grid.impl.grid.lats[l2[1]-1,l2[0]-1])
                       
                blons = np.array(blons)#.ravel().reshape(-1,2)[:,0]
                blats =  np.array(blats)#.ravel().reshape(-1,2)[:,1] 
            #                  print bound,blons,blats
                             
                tval = tide(tmodel=self.tmodel, tpath=self.tpath, blons=blons,blats=blats)
                    
            setattr(self.tide, key, tval)        
                                               
    
     
    @staticmethod 
    def to_force(ar,**kwargs):
                
        path = kwargs.get('rpath','./') 
        
        [p,u,v] = kwargs.get('vars','[None,None,None]')                
        
        curvi = kwargs.get('curvi', False)
        
        
        
        dlat=ar.lats[1,0]-ar.lats[0,0]
        dlon=ar.lons[0,1]-ar.lons[0,0]
        dlat=dlat.data.tolist() 
        dlon=dlon.data.tolist()
        lat0=ar.lats[0,0].data.tolist()
        lon0=ar.lons[0,0].data.tolist()  
        
        nodata=-9999.000
        
        pp = ar[p].fillna(nodata)
        uu = ar[u].fillna(nodata)
        vv = ar[v].fillna(nodata)
        

        if not os.path.exists(path):
           os.makedirs(path)

           # open files
        pfid = open(path+'p.amp','w')
        ufid = open(path+'u.amu','w')
        vfid = open(path+'v.amv','w')

        fi=[pfid,ufid,vfid]
        wi=[ufid,vfid]

        # write file headers
        for f in fi:
           f.write('FileVersion      = 1.03\n')
        if curvi :
           for f in fi:
              f.write('Filetype         = meteo_on_curvilinear_grid\n')
              f.write('grid_file        = wind.grd\n')
              f.write('first_data_value = grid_ulcorner\n')
              f.write('data_row         = grid_row\n')
        else:
           for f in fi:
              f.write('Filetype         = meteo_on_equidistant_grid\n')
              f.write('n_cols           = {}\n'.format(ar[u].shape[2]))
              f.write('n_rows           = {}\n'.format(ar[u].shape[1]))
              f.write('grid_unit        = degree\n')
        # code currently assumes lon and lat are increasing
              f.write('x_llcenter       = {:g}\n'.format(lon0))
              f.write('dx               = {:g}\n'.format(dlon))
              f.write('y_llcenter       = {:g}\n'.format(lat0))
              f.write('dy               = {:g}\n'.format(dlat))

        for f in fi:
           f.write('NODATA_value     = {:.3f}\n'.format(nodata))
           f.write('n_quantity       = 1\n')

        ufid.write('quantity1        = x_wind\n')
        vfid.write('quantity1        = y_wind\n')
        pfid.write('quantity1        = air_pressure\n')

        for f in wi:
           f.write('unit1            = m s-1\n')

        pfid.write('unit1            = Pa\n')

        time0=pd.to_datetime('2000-01-01 00:00:00')
    
       # write time blocks
        indx = ar.time.values - time0.to_datetime64()
        indx = indx.astype('timedelta64[m]')/60
                 
        for it in range(indx.size): # nt + 0 hour    
          for f in fi:
             f.write('TIME = {} hours since 2000-01-01 00:00:00 +00:00\n'.format(indx[it].astype(int)))

          np.savetxt(pfid,np.flipud(pp[it,:,:]),fmt='%.3f')
          np.savetxt(ufid,np.flipud(uu[it,:,:]),fmt='%.3f')
          np.savetxt(vfid,np.flipud(vv[it,:,:]),fmt='%.3f')

         # close files
        for f in fi:
           f.close()
    
            
    def run(self,**kwargs):
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
                
        bin_path = get_value(self,kwargs,'exec', None)   
            
        ncores = get_value(self,kwargs,'ncores',1)
                
        # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./run_flow2d3d.sh {} {}'.format(ncores,bin_path)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
#        for line in iter(ex.stderr.readline,b''): print line
#        ex.stderr.close() 
        with open(calc_dir+'run.log', 'w') as f: 
          for line in iter(ex.stdout.readline,b''): 
            f.write(line)   
            sys.stdout.write(line)
            sys.stdout.flush()  
        ex.stdout.close()            
                                  
    
    def pickle(self,**kwargs):
        
         path = get_value(self,kwargs,'rpath','./')  
        
         with open(path+self.tag+'.pkl', 'w') as f:
               pickle.dump(self.__dict__,f)
        
    def save(self,**kwargs):
               
         path = get_value(self,kwargs,'rpath','./')
        
#         lista = [key for key, value in self.__dict__.items() if key not in ['meteo','dem','grid']]
#         dic = {k: self.__dict__.get(k, None) for k in lista}
#         dic['model'] = self.__class__.__name__
#         dic = dict(dic, **{k: self.__dict__.get(k, None).impl.__class__.__name__ for k in ['meteo','dem']})
         with open(path+'info.pkl', 'w') as f:
               pickle.dump(self.model,f)
               
         z=self.model.copy()
         for attr, value in z.iteritems():
             if isinstance(value, datetime.datetime) : z[attr]=z[attr].isoformat()
         json.dump(z,open(path+'model.txt','w'))      
    
    @staticmethod
    def read(filename,**kwargs):
        
        filename = kwargs.get('filename', './')               
        with open(filename, 'r') as f:
              info=pickle.load(f)
        m = model(**info)
        return m
        
    def output(self,**kwargs):      
        
        path = get_value(self,kwargs,'rpath','./') 
        slevel = get_value(self,kwargs,'slevel',0.) 
        flag = get_value(self,kwargs,'update',None)
        
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        #save mdf 
        self.mdf.to_csv(path+self.tag+'.mdf',sep='=')
        
        if flag is not None:
            # check if files exist
            check=[os.path.exists(self.rpath+'{}.grd'.format(self.tag))]   
            if (np.any(check)==False) or ('dem' in flag ) :
            #save grid
                self.grid.to_file(filename = path+self.tag+'.grd')
            else:
                sys.stdout.write('skipping grid file ..\n')
        else:
            self.grid.impl.to_file(filename = path+self.tag+'.grd')
                
        
        #save dem
        try :
            bat = -self.dem.impl.fval.dem.values.astype(float) #reverse for the hydro run
       #     mask = bat==999999
        except AttributeError:    
            bat = -self.dem.impl.ival.astype(float) #reverse for the hydro run
        
        mask = ~np.isnan(bat) # mask out potential nan points
        mask[mask] = np.less(bat[mask] , 0) # get mask for dry points

        bat[mask]=np.nan #mask dry points
            
        # append the line/column of nodata 
        nodata=np.empty(self.ni)
        nodata.fill(np.nan)
        bat1=np.vstack((bat,nodata))
        nodata=np.empty((self.nj+1,1))
        nodata.fill(np.nan)
        bat2=np.hstack((bat1,nodata))

        bat2[np.isnan(bat2)] = -999.
            
        # Write bathymetry file    
        if flag is not None:
            # check if files exist
            check=[os.path.exists(self.rpath+'{}.dep'.format(self.tag))]   
            if (np.any(check)==False) or ('dem' in flag) :
                 np.savetxt(path+self.tag+'.dep',bat2)
            else:
                sys.stdout.write('skipping dem file ..\n')
        else:
            np.savetxt(path+self.tag+'.dep',bat2)
                 
        
        #save meteo
        if self.atm:
           try:
              self.to_force(self.meteo.impl.uvp,vars=['msl','u10','v10'],rpath=path,**kwargs)
           except AttributeError as e:
              print e 
              pass

             
        if self.tide :  
            
        #save bnd
            with open(path+self.tag+'.bnd', 'w') as f:
                
                dd = OrderedDict([('North',self.bound.North),('South',self.bound.South),('West',self.bound.West),('East',self.bound.East)])
            
            #    for key,val in self.bound.__dict__.iteritems():
                for i, (key, val) in enumerate(dd.iteritems()): # to match deltares 
                    
                    idx=1
                    for k1,k2 in val:           
                        bname=key+str(idx)
                        f.write('{0:<10s}{1:>12s}{2:>2s}{3:>6d}{4:>6d}{5:>6d}{6:>6d}   0.0000000e+00 {7:<s}{8:<g}A {9:<s}{10:<g}B\n'.format(bname,nm[0],nm[1],k1[0]+1,k1[1]+1,k2[0]+1,k2[1]+1,key,idx,key,idx)) # fortran index ??
                        idx+=1
        
        #save bca
            with open(path+self.tag+'.bca', 'w') as f:
            
                 dd = OrderedDict([('North',self.tide.North),('South',self.tide.South),('West',self.tide.West),('East',self.tide.East)])
            
            #     for key,val in self.tide.__dict__.iteritems():
                 for i, (key, val) in enumerate(dd.iteritems()): # to match deltares 
                     
                     idx=1
                     if val: 
                        l = np.arange(val.impl.ampl.shape[0])+idx
                        nl = [x for pair in zip(l,l) for x in pair]
                        sl = val.impl.ampl.shape[0]*le
                        for t1,t2,amp,phase in zip(np.transpose(nl),np.transpose(sl),val.impl.ampl,val.impl.phase):
                             f.write('{}{}{}\n'.format(key,t1,t2))
                             for a,b,c in zip(val.impl.constituents,amp.flatten(),phase.flatten()):
                                 f.write('{0:<3s}        {1:<.7e}   {2:<.7e}\n'.format(a,b,c))
                        
        #save obs
        if flag is not None:
        
            check=[os.path.exists(self.rpath+'{}.enc'.format(self.tag))]   
            if (np.any(check)==False) or ('model' in flag) :
            #save enc
            #write enc out
                with open(path+self.tag+'.enc','w') as f:
                    f.write('{:>5}{:>5}\n'.format(self.ni+1,1))  # add one like ddb
                    f.write('{:>5}{:>5}\n'.format(self.ni+1,self.nj+1))
                    f.write('{:>5}{:>5}\n'.format(1,self.nj+1))
                    f.write('{:>5}{:>5}\n'.format(1,1))
                    f.write('{:>5}{:>5}\n'.format(self.ni+1,1))
                
        else:
            
            #write enc out
            with open(path+self.tag+'.enc','w') as f:
                f.write('{:>5}{:>5}\n'.format(self.ni+1,1))  # add one like ddb
                f.write('{:>5}{:>5}\n'.format(self.ni+1,self.nj+1))
                f.write('{:>5}{:>5}\n'.format(1,self.nj+1))
                f.write('{:>5}{:>5}\n'.format(1,1))
                f.write('{:>5}{:>5}\n'.format(self.ni+1,1))
            
        
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
        
        try:
            bin_path = os.environ['D3D']
        except:
            bin_path = None
                
        bin_path = get_value(self,kwargs,'exec', bin_path) 
        
        if bin_path is None:
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('D3D executable path (exec) not given\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
              
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)
                
        if not os.path.exists( calc_dir+'config_d_hydro.xml') :
            
          # edit and save config file
          copy2(DATA_PATH + 'config_d_hydro.xml',calc_dir+'config_d_hydro.xml')          

        xml=md.parse(calc_dir+'config_d_hydro.xml')

        xml.getElementsByTagName('mdfFile')[0].firstChild.replaceWholeText(self.tag+'.mdf')
 
    
        with open(calc_dir+'config_d_hydro.xml','w') as f:
            xml.writexml(f)

        if not os.path.exists(calc_dir+'run_flow2d3d.sh') :

          copy2(DATA_PATH + 'run_flow2d3d.sh',calc_dir+'run_flow2d3d.sh')
        
          #make the script executable
          execf = calc_dir+'run_flow2d3d.sh'
          mode = os.stat(execf).st_mode
          mode |= (mode & 0o444) >> 2    # copy R bits to X
          os.chmod(execf, mode)
                
        
