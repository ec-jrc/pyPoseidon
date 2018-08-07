import os
import datetime
import numpy as np
import pickle
import xml.dom.minidom as md
from shutil import copy2
import subprocess
import sys
import pkg_resources
import json
from collections import OrderedDict
import pandas as pd
import glob
from shutil import copyfile
import xarray as xr

#local modules
from bnd import *
import pyPoseidon
import pyPoseidon.grid as pgrid
import pyPoseidon.meteo as pmeteo
import pyPoseidon.dem as pdem
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
        
    def config(self,**kwargs):
        self.impl.config(**kwargs)    

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

    @staticmethod               
    def read_model(filename,**kwargs):
                
        end = filename.split('.')[-1]

        if end in ['txt','json'] :
            info = pd.read_json(filename,lines=True).T.to_dict()[0]
        elif end in ['pkl']:
            info = pd.read_pickle(filename)
        
        return model(**info)
    
    @staticmethod           
    def read_folder(rfolder,**kwargs):
                
        s = glob.glob(rfolder + '/param*')
        d = glob.glob(rfolder + '/*.mdf')
                
        if s : 
            hfile = rfolder + '/hgrid.gr3' # Grid
            b = schism()
            params = pd.read_csv(s[0],engine='python',comment='!', header=None, delimiter=' =', names=['attrs','vals'])

            fix_list = [k for k in params.attrs if '=' in k ]

            for k in fix_list:
                try:
                    name, value = params.loc[params.attrs == k,'attrs'].values[0].split('=')
                    params.loc[params.attrs == k,'attrs'] = name
                    params.loc[params.attrs == name,'vals'] = value
                except:
                    pass
                
            b.params = params.set_index('attrs')
            b.grid = pgrid.grid(type='tri2d',grid_file=hfile)
                                
        elif d :
            gfile = glob.glob(rfolder + '/*.grd') # Grid
            dfile = glob.glob(rfolder + '/*.dep') # bathymetry
            
            b = d3d()
            b.mdf = pd.read_csv(d[0],sep='=')
            b.mdf = b.mdf.set_index(b.mdf.columns[0]) # set index
            b.grid = pgrid.grid('r2d',grid_file=gfile[0])
            #bath
            rdem = np.loadtxt(dfile[0])
            zx=b.grid.impl.Dataset.lons.values
            zy=b.grid.impl.Dataset.lats.values
            
            b.dem = xr.Dataset({'ival': (['ilat', 'ilon'],  rdem[:-1,:-1]), 
                             'ilons': (['k', 'l'], zx),   
                             'ilats': (['k', 'l'], zy)}, 
                             coords={'ilon': ('ilon', zx[0,:]),   
                                     'ilat': ('ilat', zy[:,0])})         
                                                            
        else:
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('No mdf nor param files present. Abort\n')
            sys.stdout.flush()
            sys.exit(1)
            #--------------------------------------------------------------------- 
        
        return b   
        
        
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
            self.time_frame = self.end_date - self.start_date
        
        if not hasattr(self, 'date'): self.date = self.start_date
        
        if self.end_date == None : 
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('model not set properly, No end_date\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
            
        
        self.tag = kwargs.get('tag', 'd3d')
        self.resolution = kwargs.get('resolution', None)
        self.ft1 = kwargs.get('ft1', 0)
        self.ft2 = kwargs.get('ft2', 11)
        self.dft = kwargs.get('dft', 1)       
        self.tide = kwargs.get('tide', False)
        self.atm = kwargs.get('atm', True)
        self.ofilename = kwargs.get('ofilename', None)

        self.solver = self.__class__.__name__    
        
        self.epath = kwargs.get('epath', None)
        
        self.Tstart = self.start_date.hour*60     
        self.Tstop = self.Tstart + int(pd.to_timedelta(self.time_frame).total_seconds()/60)

        self.step = get_value(self,kwargs,'step',0)
        self.rstep = get_value(self,kwargs,'rstep',0)
        self.dt = get_value(self,kwargs,'Dt',1)
                                               
        for attr, value in kwargs.iteritems():
                if not hasattr(self, attr): setattr(self, attr, value)
        
                          
    def set(self,**kwargs):

        gx = get_value(self,kwargs,'x',None)
        gy = get_value(self,kwargs,'y',None)    
              
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
          
        self.grid=pgrid.grid(type='r2d',x=gx, y=gy)
        
        # get bathymetry
        self.bath(**kwargs)

        # get boundaries
        self.bc()
                
        #get meteo
        if self.atm :  self.force(**kwargs)
        
        #get tide
        if self.tide : self.tidebc()
        
        self.config(**kwargs)
        
        #mdf
    def config(self,**kwargs):  
        
        mdf_file = kwargs.get('config_file', None)  

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
        if self.ofilename :
            self.mdf.loc[self.mdf.index.str.contains('Filsta')]='#{}#'.format(self.tag+'.obs') 
        else:
            self.mdf.loc[self.mdf.index.str.contains('Filsta')]='##'
   
        # adjust ni,nj
        self.mdf.loc[self.mdf.index.str.contains('MNKmax')]='{} {} {}'.format(self.ni+1,self.nj+1,1)  # add one like ddb
  
        # adjust iteration date
        self.mdf.loc[self.mdf.index.str.contains('Itdate')]='#{}#'.format(self.date.strftime(format='%Y-%m-%d'))
  
        #set time unit
        self.mdf.loc[self.mdf.index.str.contains('Tunit')]='#M#'

        #adjust iteration start
        self.mdf.loc[self.mdf.index.str.contains('Tstart')]=self.Tstart
  
        #adjust iteration stop
        self.mdf.loc[self.mdf.index.str.contains('Tstop')]=self.Tstop
  
        #adjust time step
        self.mdf.loc[self.mdf.index.str.contains('Dt')]=[self.dt]
  
        #adjust time for output
        self.mdf.loc[self.mdf.index.str.contains('Flmap')]='{} {} {}'.format(self.Tstart,self.step,self.Tstop)
        self.mdf.loc[self.mdf.index.str.contains('Flhis')]='{} {} {}'.format(self.Tstart,self.dt,self.Tstop)
        self.mdf.loc[self.mdf.index.str.contains('Flpp')]='0 0 0'
        self.mdf.loc[self.mdf.index.str.contains('Flrst')]=self.rstep
  
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

        other = kwargs.get('config', None) 
        if other:
        # Check for any other mdf variable in input
            for key,val in other.iteritems():
                if key in mdfidx: 
                    self.mdf.loc[self.mdf.index.str.contains(key)] = val
                else:
                    self.mdf.loc[key] = val
        
        
        output = kwargs.get('output', False)
        
        if output: 
            #save mdf 
            path = get_value(self,kwargs,'rpath','./') 
            self.mdf.to_csv(path+self.tag+'.mdf',sep='=')
            

        #meteo
    def force(self,**kwargs):
        z = self.__dict__.copy()
                
        mpaths =  get_value(self,kwargs,'mpaths',None)        
        
        z.update({'mpaths':mpaths})

        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:     
            check=[os.path.exists(z['rpath']+f) for f in ['u.amu','v.amv','p.amp']]   
            if (np.any(check)==False) :              
                self.meteo = pmeteo.meteo(**z)
            elif 'meteo' in flag : 
                self.meteo = pmeteo.meteo(**z)        
            else:
                sys.stdout.write('skipping meteo files ..\n')
        else:
            self.meteo = pmeteo.meteo(**z)
        
        #dem
    def bath(self,**kwargs):
        z = self.__dict__.copy()        
        
        z['grid_x'] = self.grid.impl.Dataset.lons.values
        z['grid_y'] = self.grid.impl.Dataset.lats.values
        
        dpath =  get_value(self,kwargs,'dpath',None)        
        
        z.update({'dpath':dpath})
        
        if 'dem_file' in kwargs.keys():
            
            rdem = np.loadtxt(kwargs['dem_file'])
            
            self.dem = xr.Dataset({'ival': (['ilat', 'ilon'],  rdem[:-1,:-1]), 
                             'ilons': (['k', 'l'], z['grid_x']),   
                             'ilats': (['k', 'l'], z['grid_y'])}, 
                             coords={'ilon': ('ilon', z['grid_x'][0,:]),   
                                     'ilat': ('ilat', z['grid_y'][:,0])})         
            
        
        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:
            check=[os.path.exists(z['rpath']+f) for f in ['{}.dep'.format(z['tag'])]]   
            if (np.any(check)==False) :
                self.dem = pdem.dem(**z)  
            elif 'dem' in flag :
                self.dem = pdem.dem(**z)
            else:
                sys.stdout.write('reading local dem file ..\n')
                dem_file = z['rpath']+self.tag+'.dep'
                rdem = np.loadtxt(dem_file)
                zx=self.grid.impl.Dataset.lons.values
                zy=self.grid.impl.Dataset.lats.values
            
                self.dem = xr.Dataset({'ival': (['ilat', 'ilon'],  rdem[:-1,:-1]), 
                                 'ilons': (['k', 'l'], zx),   
                                 'ilats': (['k', 'l'], zy)}, 
                                 coords={'ilon': ('ilon', zx[0,:]),   
                                         'ilat': ('ilat', zy[:,0])})         
                
        else:
            self.dem = pdem.dem(**z)
                
        
    def bc(self,**kwargs):
        #define boundaries
        z = self.__dict__.copy()        
        
        z['lons'] = self.grid.impl.Dataset.lons[0,:]
        z['lats'] = self.grid.impl.Dataset.lats[:,0]
        
        try:
            ba = -self.dem.impl.Dataset.ival.astype(np.float)
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
                    blons.append(self.grid.impl.Dataset.lons[l1[1]-1,l1[0]-1])   
                    blats.append(self.grid.impl.Dataset.lats[l1[1]-1,l1[0]-1])
                    blons.append(self.grid.impl.Dataset.lons[l2[1]-1,l2[0]-1])   
                    blats.append(self.grid.impl.Dataset.lats[l2[1]-1,l2[0]-1])
                       
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
                
        bin_path = get_value(self,kwargs,'epath', None)   
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)
        
        conda_bin = get_value(self,kwargs,'conda_bin', None)
        
        argfile = get_value(self,kwargs,'argfile',self.tag+'_hydro.xml')
                        
        if conda_env is None:        
            # note that cwd is the folder where the executable is
            ex=subprocess.Popen(args=['./run_flow2d3d.sh {} {} {}'.format(argfile,ncores,bin_path)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
        else:
            ex=subprocess.Popen(args=['./run_flow2d3d.sh {} {} {} {} {}'.format(argfile,ncores,bin_path,conda_bin,conda_env)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
        
        
        with open(calc_dir+self.tag+'_run.log', 'w') as f: #save output
            
            for line in iter(ex.stdout.readline,b''): 
                f.write(line)   
                sys.stdout.write(line)
                sys.stdout.flush()  
            
            for line in iter(ex.stderr.readline,b''):
                sys.stdout.write(line)
                sys.stdout.flush()  
                tempfiles = glob.glob(calc_dir+'/tri-diag.'+ self.tag+'-*')
                try:
                    biggest = max(tempfiles, key=(lambda tf: os.path.getsize(tf)))
                    with open(biggest, "r") as f1:
                        for line in f1:
                            f.write(line)
                except:
                    pass
                            
      #cleanup  
        tempfiles = glob.glob(calc_dir+'/tri-diag.'+ self.tag+'-*') + glob.glob(calc_dir+'/TMP_*')
                                
        for filename in tempfiles:
            try:
                os.remove(filename)
            except OSError:
                pass

        ex.stdout.close()  
        ex.stderr.close() 

    
    def pickle(self,**kwargs):
        
         path = get_value(self,kwargs,'rpath','./')  
        
         with open(path+self.tag+'.pkl', 'w') as f:
               pickle.dump(self.__dict__,f)
        
    def save(self,**kwargs):
               
         path = get_value(self,kwargs,'rpath','./')
        
         lista = [key for key, value in self.__dict__.items() if key not in ['meteo','dem','grid']]
         dic = {k: self.__dict__.get(k, None) for k in lista}

         grid=self.__dict__.get('grid', None)
         if isinstance(grid,np.str):
             dic.update({'grid':grid})
         elif isinstance(grid,pgrid.grid):
             dic.update({'grid':grid.impl.__class__.__name__})
         
         dem=self.__dict__.get('dem', None)
         if isinstance(dem,np.str):
             dic.update({'dem':dem})
         elif isinstance(dem,pdem.dem):
             dic.update({'dem':dem.impl.__class__.__name__})

         meteo=self.__dict__.get('meteo', None)
         if isinstance(meteo,np.str):
             dic.update({'meteo':meteo})
         elif isinstance(meteo,pmeteo.meteo):
             dic.update({'meteo':meteo.impl.__class__.__name__})

         dic['version']=pyPoseidon.__version__
         
         with open(path+self.tag+'_info.pkl', 'w') as f:
               pickle.dump(dic,f)
                        
         for attr, value in dic.iteritems():
             if isinstance(value, datetime.datetime) : dic[attr]=dic[attr].isoformat()
             if isinstance(value, pd.Timedelta) : dic[attr]=dic[attr].isoformat()          
             if isinstance(value, pd.DataFrame) : dic[attr]=dic[attr].to_dict()
         json.dump(dic,open(path+self.tag+'_model.json','w'))      
    
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
            if (np.any(check)==False) or ('grid' in flag ) :
            #save grid
                self.grid.impl.to_file(filename = path+self.tag+'.grd')
            else:
                sys.stdout.write('skipping grid file ..\n')
        else:
            self.grid.impl.to_file(filename = path+self.tag+'.grd')
                
        
        #save dem
        try:
            try :
                bat = -self.dem.impl.Dataset.fval.values.astype(float) #reverse for the hydro run
       #     mask = bat==999999
            except AttributeError:    
                bat = -self.dem.impl.Dataset.ival.values.astype(float) #reverse for the hydro run
        
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
                        
        except AttributeError:
            sys.stdout.write('No dem file ..')
            
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
             # print e 
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
        
        ofilename = get_value(self,kwargs,'ofilename',None) 

        if ofilename:
        
            obs_points = pd.read_csv(ofilename,delimiter='\t',header=None,names=['index','Name','lat','lon'])
            obs_points = obs_points.set_index('index',drop=True).reset_index(drop=True) #reset index if any
        
            b=np.ma.masked_array(bat,np.isnan(bat)) # mask land
        
            i_indx, j_indx = self.vpoints(self.grid.impl.Dataset,obs_points,b,**kwargs)
                    
            obs_points['i']=i_indx
            obs_points['j']=j_indx
        
            #drop NaN points
            obs = obs_points.dropna().copy()
        
            obs = obs.reset_index(drop=True) #reset index
        
            obs['i']=obs['i'].values.astype(int)
            obs['j']=obs['j'].values.astype(int)
            obs['new_lat']=self.grid.impl.Dataset.y[obs.i.values].values #Valid point
            obs['new_lon']=self.grid.impl.Dataset.x[obs.j.values].values 
        
            self.obs = obs #store it
          
            obs.Name = obs.Name.str.strip().apply(lambda name:name.replace(' ', '')) #Remove spaces to write to file
            sort = sorted(obs.Name.values,key=len) # sort the names to get the biggest word
            wsize = len(sort[-1])# size of bigget word in order to align below
        
            if flag is not None:
        
                check=[os.path.exists(self.rpath+'{}.obs'.format(self.tag))]   
                if (np.any(check)==False) or ('model' in flag) :
        
                    # Add one in the indices due to python/fortran convention
                    with open(self.rpath+'{}.obs'.format(self.tag),'w') as f: 
                        for l in range(obs.shape[0]): 
                            f.write('{0:<{3}}{1:>{3}}{2:>{3}}\n'.format(obs.Name[l],obs.j[l]+1,obs.i[l]+1,wsize))

         
            else:
                # Add one in the indices due to python/fortran convention
                with open(self.rpath+'{}.obs'.format(self.tag),'w') as f:
                    for l in range(obs.shape[0]): 
                        f.write('{0:<{3}}{1:>{3}}{2:>{3}}\n'.format(obs.Name[l],obs.j[l]+1,obs.i[l]+1,wsize))
              
        
        #save enc file
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
                        
        bin_path = get_value(self,kwargs,'epath', None) 
        
        if bin_path is None:
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('D3D executable path (epath) not given\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
              
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)
                        
        if not os.path.exists( calc_dir+self.tag+'_hydro.xml') :
            
          # edit and save config file
          copy2(DATA_PATH + 'config_d_hydro.xml',calc_dir+self.tag+'_hydro.xml')          

        xml=md.parse(calc_dir+self.tag+'_hydro.xml')

        xml.getElementsByTagName('mdfFile')[0].firstChild.replaceWholeText(self.tag+'.mdf')
 
    
        with open(calc_dir+self.tag+'_hydro.xml','w') as f:
            xml.writexml(f)

        if not os.path.exists(calc_dir+'run_flow2d3d.sh') :

          copy2(DATA_PATH + 'run_flow2d3d.sh',calc_dir+'run_flow2d3d.sh')
        
          #make the script executable
          execf = calc_dir+'run_flow2d3d.sh'
          mode = os.stat(execf).st_mode
          mode |= (mode & 0o444) >> 2    # copy R bits to X
          os.chmod(execf, mode)
                
    @staticmethod
    def vpoints(grid,obs_points,bat,**kwargs):
        
        idx=[]
        jdx=[]
        for m in range(obs_points.shape[0]):
            lat, lon = obs_points.loc[m,['lat','lon']]
            nearest = grid.sel(x=[lon],y=[lat], method='nearest')
            j = np.abs(grid.x.values-nearest.x.values).argmin()
            i = np.abs(grid.y.values-nearest.y.values).argmin()
            if bat[i,j] :
                idx.append(i) 
                jdx.append(j)
            else:
                bnear=bat[i-5:i+6,j-5:j+6] # near by grid nodes
        
                rlon = grid.lons[i-5:i+6,j-5:j+6]-lon
                rlat = grid.lats[i-5:i+6,j-5:j+6]-lat
                rad = np.sqrt(rlon**2+rlat**2) # radial distance from the obs point
        
                rmask = rad.values[bnear.mask==False] #mask the distance array with the valid mask from dem
                
                rmask.sort() # sort to start close and move further away
                if rmask.size > 0 :
                    
                    for r in rmask: # Find the closest valid point
                        [[k,l]] = np.argwhere(rad.values==r)
                        if bnear[k-1:k+1,l-1:l+1].mask.sum() == 0:
                            break# The point is valid point

                    xv = rad[k,l].x.values #lat, lon of valid point
                    yv = rad[k,l].y.values
        
                    #final i,j
                    j = np.abs(grid.x.values-xv).argmin()
                    i = np.abs(grid.y.values-yv).argmin()

                    idx.append(i) 
                    jdx.append(j) 
          
                else:
            
                    idx.append(np.nan)
                    jdx.append(np.nan)

        return idx,jdx
        
class schism(model):
    
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
        elif 'end_date' in kwargs:
            end_date = kwargs.get('end_date', None)
            self.end_date = pd.to_datetime(end_date)
            self.time_frame = self.end_date - self.start_date

        if not hasattr(self, 'date'): self.date = self.start_date
        
        if self.end_date == None : 
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('model not set properly, No end_date\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
            
        
        self.tag = kwargs.get('tag', 'schism')
        self.resolution = kwargs.get('resolution', None)
        self.ft1 = kwargs.get('ft1', 0)
        self.ft2 = kwargs.get('ft2', 11)
        self.dft = kwargs.get('dft', 1)       
        self.tide = kwargs.get('tide', False)
        self.atm = kwargs.get('atm', True)
    
        self.epath = kwargs.get('epath', None)
    
        self.Tstart = self.start_date.hour*60     
        self.Tstop = self.Tstart + int((self.end_date - self.start_date).total_seconds()/60)
        
        self.step = get_value(self,kwargs,'step',0)
        self.rstep = get_value(self,kwargs,'rstep',0)
        self.dt = get_value(self,kwargs,'dt',400)
        
        self.solver = self.__class__.__name__    
        
                                                   
        for attr, value in kwargs.iteritems():
            if not hasattr(self, attr): setattr(self, attr, value)  
         
                          
    def set(self,**kwargs):


        self.hgrid = kwargs.get('grid_file',None)
                                         
        # Grid         
        self.grid=pgrid.grid(type='tri2d',**kwargs)
                 
        # set lat/lon from file
        if self.hgrid:
            self.minlon = self.grid.impl.Dataset.x.dropna('id').values.min()
            self.maxlon = self.grid.impl.Dataset.x.dropna('id').values.max()
            self.minlat = self.grid.impl.Dataset.y.dropna('id').values.min()
            self.maxlat = self.grid.impl.Dataset.y.dropna('id').values.max()
                                     
        # get bathymetry
        self.bath(**kwargs)

        # get boundaries
        # self.bc()
                
        #get meteo
        if self.atm :  self.force(**kwargs)
        
        #get tide
        if self.tide : self.tidebc()
        
        
        self.config(**kwargs)
    
        #param
    def config(self,**kwargs): 
        
        dic = get_value(self,kwargs,'configx    ',None)
        param_file = get_value(self,kwargs,'config_file',None)
           
        if param_file :
            params = pd.read_csv(param_file,engine='python',comment='!', header=None, delimiter=' =', names=['attrs','vals'])

            fix_list = [k for k in params.attrs if '=' in k ]

            for k in fix_list:
                try:
                    name, value = params.loc[params.attrs == k,'attrs'].values[0].split('=')
                    params.loc[params.attrs == k,'attrs'] = name
                    params.loc[params.attrs == name,'vals'] = value
                except:
                    pass
                
            params = params.set_index('attrs')
        else:
            with open(DATA_PATH + 'dparams.pkl', 'r') as f:
                              params=pickle.load(f)
        
        if dic :
            for key,val in dic.iteritems():
                params.loc[key] = val
            
        #update 
        params.loc['start_hour'] = self.start_date.hour
        params.loc['start_day'] = self.start_date.day
        params.loc['start_month'] = self.start_date.month
        params.loc['start_year'] = self.start_date.year
        params.loc['rnday'] = (self.Tstop-self.Tstart)/(60*24.)
        params.loc['dt'] = self.dt
            
        
        self.params = params
        
        output = kwargs.get('output', False)
                
        if output: 
            #save params 
            path = get_value(self,kwargs,'rpath','./') 
            self.params.to_csv(path + 'param.in', header=None, sep='=')  #save to file

        #meteo
    def force(self,**kwargs):
        z = self.__dict__.copy()
        
        mpaths =  get_value(self,kwargs,'mpaths',None)        
        
        z.update({'mpaths':mpaths})
        
        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:     
            check=[os.path.exists(z['rpath']+'sflux/sflux_air_1.001.nc')] # for f in ['u.amu','v.amv','p.amp']]   
            if (np.any(check)==False) :              
                self.meteo = pmeteo.meteo(**z)
            elif 'meteo' in flag : 
                self.meteo = pmeteo.meteo(**z)        
            else:
                sys.stdout.write('skipping meteo ..\n')
        else:
            self.meteo = pmeteo.meteo(**z)
        
        #dem
    def bath(self,**kwargs):
        z = self.__dict__.copy()        
        
        z['grid_x'] = self.grid.impl.Dataset.x.dropna('id').values
        z['grid_y'] = self.grid.impl.Dataset.y.dropna('id').values
        
        dpath =  get_value(self,kwargs,'dpath',None)        
        
        z.update({'dpath':dpath})
                
        flag = get_value(self,kwargs,'update',None)
        # check if files exist
        if flag is not None:
            if 'dem' in flag :
                self.dem = pdem.dem(**z)
            else:
                sys.stdout.write('dem from grid file\n')
        else:
            self.dem = pdem.dem(**z)

    def output(self,**kwargs):      
        
        path = get_value(self,kwargs,'rpath','./') 
        flag = get_value(self,kwargs,'update',None)
        
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        # save sflux_inputs.txt
        if not os.path.exists(path+'sflux'):
            os.makedirs(path+'sflux')
        
        with open(path+'sflux/sflux_inputs.txt', 'w') as f:
            f.write('&sflux_inputs\n')
            f.write('/ \n\n')
            
        # save bctides.in
        nobs = [key for key in self.grid.impl.Dataset.variables.keys() if 'open' in key]
        
        with open(path + 'bctides.in', 'w') as f:
            f.write('Header\n')
            f.write('{} {}\n'.format(0, 40.)) #  ntip tip_dp
            f.write('{}\n'.format(0)) #nbfr
            f.write('{}\n'.format(len(nobs))) #number of open boundaries
            for i in range(len(nobs)):
                nnodes = self.grid.impl.Dataset.to_dataframe().loc[:,[nobs[i]]].dropna().size
                f.write('{} {} {} {} {}\n'.format(nnodes,2,0,0,0)) # number of nodes on the open boundary segment j (corresponding to hgrid.gr3), B.C. flags for elevation, velocity, temperature, and salinity
                f.write('{}\n'.format(0)) # ethconst !constant elevation value for this segment    
       
       #save vgrid.in
        with open(path + 'vgrid.in', 'w') as f:
            f.write('{}\n'.format(2)) #ivcor (1: LSC2; 2: SZ)
            f.write('{} {} {}\n'.format(2,1,1.e6)) #nvrt(=Nz); kz (# of Z-levels); hs (transition depth between S and Z)
            f.write('Z levels\n') # Z levels !Z-levels in the lower portion
            f.write('{} {}\n'.format(1,-1.e6)) #!level index, z-coordinates, z-coordinate of the last Z-level must match -hs 
            f.write('S levels\n') # S-levels below 
            f.write('{} {} {}\n'.format(40.,1.,1.e-4))  #constants used in S-transformation: h_c, theta_b, theta_f
            f.write('{} {}\n'.format(1,-1.)) #first S-level (sigma-coordinate must be -1)
            f.write('{} {}\n'.format(2,0.)) # levels index, sigma-coordinate, last sigma-coordinate must be 0

       # save params.in
        
        self.params.to_csv(path + 'param.in', header=None, sep='=')  #save to file
             
        
        #save hgrid.gr3
        try:
            
            bat = -self.dem.impl.Dataset.ival.values.astype(float) #minus for the hydro run
                                    
            self.grid.impl.Dataset.z.loc[:bat.size] = bat
                                    
            self.grid.impl.to_file(filename= path+'hgrid.gr3')
            
        except AttributeError as e:
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('Keeping bathymetry from hgrid.gr3 ..\n')
            sys.stdout.flush()
            
            
        # save grid files 
        
        copyfile(self.hgrid, path+'hgrid.gr3') #copy original grid file
        
        copyfile(path+'hgrid.gr3', path+'hgrid.ll')    
                 
                 
        # manning file
        manfile=path+'manning.gr3'
        
        nn = self.grid.impl.Dataset.x[np.isfinite(self.grid.impl.Dataset.x.values)].size
        n3e = self.grid.impl.Dataset.a.size
        
        with open(manfile,'w') as f:
            f.write('\t 0 \n')
            f.write('\t {} {}\n'.format(n3e,nn))
        
        mn = xr.Dataset({'man': ('id', np.ones(nn)*.12), 'id': np.arange(1,nn+1)})
                    
        self.grid.impl.Dataset=xr.merge([self.grid.impl.Dataset, mn])
        
        man = self.grid.impl.Dataset.get(['x','y','man']).to_dataframe().dropna()
        
        man.to_csv(manfile,index=True, sep='\t', header=None,mode='a', float_format='%.10f',columns=['x','y','man'] )
        
        e = self.grid.impl.Dataset.get(['nv','a','b', 'c']).to_dataframe()
        
        e.to_csv(manfile,index=True, sep='\t', header=None, mode='a', columns=['nv','a','b','c']) 
        
        
        # windrot_geo2proj
        
        windfile=path+'windrot_geo2proj.gr3'
        
        with open(windfile,'w') as f:
            f.write('\t 0 \n')
            f.write('\t {} {}\n'.format(n3e,nn))
            
        wdr = xr.Dataset({'windrot': ('id', np.ones(nn)*0.000001), 'id': np.arange(1,nn+1)})
        
        self.grid.impl.Dataset=xr.merge([self.grid.impl.Dataset, wdr])
               
        wd = self.grid.impl.Dataset.get(['x','y','windrot']).to_dataframe().dropna()
        
        wd.to_csv(windfile,index=True, sep='\t', header=None,mode='a', float_format='%.10f',columns=['x','y','windrot'] )
        
        e.to_csv(windfile,index=True, sep='\t', header=None, mode='a', columns=['nv','a','b','c'])
        
        
        #save meteo
        if self.atm:
           try:
              self.to_force(self.meteo.impl.uvp,vars=['msl','u10','v10'],rpath=path,**kwargs)
           except AttributeError as e:
              print e 
              pass

             
                                
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
                        
        bin_path = get_value(self,kwargs,'epath', None) 
        
        if bin_path is None:
            #--------------------------------------------------------------------- 
            sys.stdout.flush()
            sys.stdout.write('\n')
            sys.stdout.write('Schism executable path (epath) not given\n')
            sys.stdout.flush()
            #--------------------------------------------------------------------- 
              
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)
                        
            
        with open(calc_dir + 'launchSchism.sh', 'w') as f:
            if conda_env :
                f.write('source activate {}\n'.format(conda_env))
            f.write('exec={}\n'.format(bin_path))
            f.write('mkdir outputs\n')
            f.write('mpirun -N {} $exec\n'.format(ncores))   
                 
          #make the script executable
        execf = calc_dir+'launchSchism.sh'
        mode = os.stat(execf).st_mode
        mode |= (mode & 0o444) >> 2    # copy R bits to X
        os.chmod(execf, mode)

    @staticmethod 
    def to_force(ar,**kwargs):
                
        path = kwargs.get('rpath','./') 
        
        [p,u,v] = kwargs.get('vars','[None,None,None]')                
            
        xx, yy = np.meshgrid(ar.longitude.data, ar.latitude.data) 
        
        zero = np.zeros(ar[p].data.shape)
                       
        udate = pd.to_datetime(ar.time[0].data).strftime('%Y-%m-%d')
        
        bdate = pd.to_datetime(ar.time[0].data).strftime('%Y %m %d %H').split(' ')
        
        tlist = (ar.time.data - pd.to_datetime([udate]).values).astype('timedelta64[s]')/3600.
        
        tlist = tlist.astype(float)/24.
        
        
        bdate = [int(q) for q in bdate]
        
        sout= xr.Dataset({'prmsl':(['time', 'nx_grid', 'ny_grid'], ar[p].data),
                          'uwind':(['time','nx_grid','ny_grid'], ar[u].data),
                          'vwind':(['time','nx_grid','ny_grid'], ar[v].data),
                          'spfh':(['time','nx_grid','ny_grid'], zero),
                          'stmp':(['time','nx_grid','ny_grid'], zero),
                          'lon':(['nx_grid','ny_grid'], xx),
                          'lat':(['nx_grid','ny_grid'], yy)},
                     coords={'time':tlist})
                     
        sout.attrs={'description' : 'Schism forsing',
            'history' :'JRC Ispra European Commission',
            'source' : 'netCDF4 python module'}
            
        sout.time.attrs={   'long_name':      'Time',
                            'standard_name':  'time',
                            'base_date':      bdate,
                            'units':          udate }
                            
        sout.lat.attrs={'units': 'degrees_north',
                       'long_name': 'Latitude',
                       'standard_name':'latitude'}
                       
        sout.prmsl.attrs={'units': 'Pa',
                       'long_name': 'Pressure reduced to MSL',
                       'standard_name':'air_pressure_at_sea_level'}
                       
        sout.uwind.attrs={'units': 'm/s',
                       'long_name': 'Surface Eastward Air Velocity',
                       'standard_name':'eastward_wind'}
                       
        sout.vwind.attrs={'units': 'm/s',
                       'long_name': 'Surface Northward Air Velocity',
                       'standard_name':'northward_wind'}
                       
        sout.spfh.attrs={'units': '1',
                       'long_name': 'Surface Specific Humidity (2m AGL)',
                       'standard_name':'specific_humidity'}               
                       
        sout.stmp.attrs={'units': 'degrees',
                       'long_name': 'Surface Temperature',
                       'standard_name':'surface temperature'}
               
        
        #check if folder sflux exists
        if not os.path.exists(path+'sflux'):
            os.makedirs(path+'sflux')
               
        sout.to_netcdf(path+'sflux/sflux_air_1.001.nc')
                                         
                
    def run(self,**kwargs):
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
                
        bin_path = get_value(self,kwargs,'epath', None)   
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)
                
            # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./launchSchism.sh'], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
            
        with open(calc_dir+'err.log', 'w') as f: 
          for line in iter(ex.stderr.readline,b''): 
            f.write(line)   
            sys.stdout.write(line)
            sys.stdout.flush()  
        ex.stderr.close()            

        with open(calc_dir+'run.log', 'w') as f: 
          for line in iter(ex.stdout.readline,b''): 
            f.write(line)   
            sys.stdout.write(line)
            sys.stdout.flush()  
        ex.stdout.close()         
        
                                  
    def save(self,**kwargs):
               
         path = get_value(self,kwargs,'rpath','./')
        
         lista = [key for key, value in self.__dict__.items() if key not in ['meteo','dem','grid']]
         dic = {k: self.__dict__.get(k, None) for k in lista}

         grid=self.__dict__.get('grid', None)
         if isinstance(grid,np.str):
             dic.update({'grid':grid})
         elif isinstance(grid,pgrid.grid):
             dic.update({'grid':grid.impl.__class__.__name__})
         
         dem=self.__dict__.get('dem', None)
         if isinstance(dem,np.str):
             dic.update({'dem':dem})
         elif isinstance(dem,pdem.dem):
             dic.update({'dem':dem.impl.__class__.__name__})

         meteo=self.__dict__.get('meteo', None)
         if isinstance(meteo,np.str):
             dic.update({'meteo':meteo})
         elif isinstance(meteo,pmeteo.meteo):
             dic.update({'meteo':meteo.impl.__class__.__name__})

         dic['version']=pyPoseidon.__version__

         with open(path+self.tag+'_info.pkl', 'w') as f:
               pickle.dump(dic,f)
               
         for attr, value in dic.iteritems():
             if isinstance(value, datetime.datetime) : dic[attr]=dic[attr].isoformat()
             if isinstance(value, pd.Timedelta) : dic[attr]=dic[attr].isoformat()          
             if isinstance(value, pd.DataFrame) : dic[attr]=dic[attr].to_dict()
         json.dump(dic,open(path+self.tag+'_model.json','w'))      
