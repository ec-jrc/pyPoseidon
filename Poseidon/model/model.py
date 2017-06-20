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
#local modules
import mdf
from grid import *
from dep import *
from Poseidon.meteo import meteo
from Poseidon.dem import dem


#retrieve the module path
DATA_PATH = pkg_resources.resource_filename('Poseidon', 'misc/')
#DATA_PATH = os.path.dirname(Poseidon.__file__)+'/misc/'    
#info_data = ('lon0','lon1','lat0','lat1','date','tag','resolution','ft1','ft2')

class model:
    impl = None
    def __init__(self, model=None, **kwargs):
        if model == 'd3d':
            self.impl = d3d(**kwargs)
        elif model == 'schism':
            self.impl = schism(**kwargs)            
        
    def create(self,**kwargs):
        self.impl.create(**kwargs)
        
    def force(self,**kwargs):
        self.impl.force(**kwargs)
                  
    def output(self,**kwargs):
        self.impl.output(**kwargs)
            
    def run(self,**kwargs):
        self.impl.run(**kwargs)
        
    def save(self,**kwargs):
        self.impl.save(**kwargs)
    
    def read(self,**kwargs):
        self.impl.read(**kwargs)

    def set(self,**kwargs):
        self.impl.set(**kwargs)
        
    def meteo(self,**kwargs):
        self.impl.meteo(**kwargs)

    def dem(self,**kwargs):
        self.impl.dem(**kwargs)
                                   
        
class d3d(model):
    
    def __init__(self,**kwargs):
        
        self.lon0 = kwargs.get('lon0', None)
        self.lon1 = kwargs.get('lon1', None)
        self.lat0 = kwargs.get('lat0', None)
        self.lat1 = kwargs.get('lat1', None)       
        self.date = kwargs.get('date', None)
        self.tag = kwargs.get('tag', None)
        self.resolution = kwargs.get('resolution', None)
        self.ft1 = kwargs.get('ft1', None)
        self.ft2 = kwargs.get('ft2', None)
        
            
    def create(self,**kwargs):

        gx = kwargs.get('x', None)
        gy = kwargs.get('y', None)    
        mdf_file = kwargs.get('mdf', None)  
        Tstart = self.date.hour+self.ft1*60     
        Tstop = self.date.hour+self.ft2*60
        step = kwargs.get('step', None)
        rstep = kwargs.get('rstep', None)
                                
        resmin=self.resolution*60
      
        # computei ni,nj / correct lat/lon

        if gx :
            
          self.x=gx
          self.y=gy
              
          nj,ni=self.x.shape
          self.lon0=self.x.min()
          self.lon1=self.x.max()
          self.lat0=self.y.min()
          self.lat1=self.y.max()

        else:

          ni=int(round((self.lon1-self.lon0)/self.resolution)) #these are cell numbers
          nj=int(round((self.lat1-self.lat0)/self.resolution))
  
          self.lon1=self.lon0+ni*self.resolution #adjust max lon to much the grid
          self.lat1=self.lat0+nj*self.resolution

       
        ni=ni+1 # transfrom to grid points
        nj=nj+1
        
        self.ni=ni
        self.nj=nj
        
        if gx is None :
        # set the grid 
          x=np.linspace(self.lon0,self.lon1,self.ni)
          y=np.linspace(self.lat0,self.lat1,self.nj)
          gx,gy=np.meshgrid(x,y)
          
        # Grid   
        self.grid=Grid()
        self.grid.x,self.grid.y=gx,gy
        self.grid.shape = gx.shape        
        
        #mdf
        if mdf_file :
            inp, order = mdf.read(mdf_file)
        else:
            inp, order = mdf.read(DATA_PATH+'default.mdf')
        
            self.mdf = Bunch({'inp':inp, 'order':order})    
            
            #define grid file
            self.mdf.inp['Filcco']=self.tag+'.grd'
  
            #define enc file
            self.mdf.inp['Filgrd']=self.tag+'.enc'
  
            #define dep file
            self.mdf.inp['Fildep']=self.tag+'.dep'
  
            #define obs file
            self.mdf.inp['Filsta']='' #b.tag+'.obs'
   
            # adjust ni,nj
            self.mdf.inp['MNKmax']=[self.ni+1,self.nj+1,1]  # add one like ddb
  
            # adjust iteration date
            self.mdf.inp['Itdate']=datetime.datetime.strftime(self.date,'%Y-%m-%d')
  
            #set time unit
            self.mdf.inp['Tunit']='M'

            #adjust iteration start
            self.mdf.inp['Tstart']=[Tstart]
  
            #adjust iteration stop
            self.mdf.inp['Tstop']=[Tstop]
  
            #adjust time step
            self.mdf.inp['Dt']=[1.]
  
            #adjust time for output
            self.mdf.inp['Flmap']=[Tstart,step,Tstop]
            self.mdf.inp['Flhis']=[Tstart,1,Tstop]
            self.mdf.inp['Flpp']=[0,0,0]
            self.mdf.inp['Flrst']=[rstep]
  
            #time interval to smooth the hydrodynamic boundary conditions
            self.mdf.inp['Tlfsmo']=[0.]


          # set tide only run
        #if 'tide' in kwargs.keys() :
        #     if (kwargs['tide']==True) & (force==False):
        #        inp['Sub1'] = ' '
        #        inp['Filbnd']=basename+'.bnd'
        #        inp['Filana']=basename+'.bca'
        #       if 'Tidfor' not in order: order.append('Tidfor')
        #       inp['Tidfor']=[['M2','S2','N2','K2'], \
        #                      ['K1','O1','P1','Q1'], \
        #                      ['-----------']]
 
          # specify ini file
        # if 'Filic' not in order: order.append('Filic')
        # inp['Filic']=basename+'.ini'

          # netCDF output
        if 'FlNcdf' not in self.mdf.order: self.mdf.order.append('FlNcdf')

        self.mdf.inp['FlNcdf'] = 'map his'

               
        #meteo
    def meteo(self,**kwargs):
        self.meteo = meteo(**kwargs)  
        
        #dem
    def dem(self,**kwargs):
        self.dem = dem(**kwargs)  

    
    def force(self,**kwargs):
    
        path = kwargs.get('rpath', './')
        curvi = kwargs.get('curvi', False)
        
        dlat=self.meteo.impl.lats[1,0]-self.meteo.impl.lats[0,0]
        dlon=self.meteo.impl.lons[0,1]-self.meteo.impl.lons[0,0]
        lat0=self.meteo.impl.lats[0,0] 
        lon0=self.meteo.impl.lons[0,0] 

        nodata=-9999.000

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
              f.write('n_cols           = {}\n'.format(self.meteo.impl.u.shape[2]))
              f.write('n_rows           = {}\n'.format(self.meteo.impl.u.shape[1]))
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

        time0=datetime.datetime.strptime('2000-01-01 00:00:00','%Y-%m-%d %H:%M:%S')

       # write time blocks
        for it in range(self.ft1,self.ft2+1): # nt + 0 hour    
          ntime=self.date+datetime.timedelta(hours=it)
          dt=(ntime-time0).total_seconds()/3600.


          for f in fi:
             f.write('TIME = {} hours since 2000-01-01 00:00:00 +00:00\n'.format(dt))

          np.savetxt(pfid,np.flipud(self.meteo.impl.p[it,:,:]/0.01),fmt='%.3f')
          np.savetxt(ufid,np.flipud(self.meteo.impl.u[it,:,:]),fmt='%.3f')
          np.savetxt(vfid,np.flipud(self.meteo.impl.v[it,:,:]),fmt='%.3f')

    # write the same values for the end time
    #dt=dt+nt*60.
    #for f in fi:
    # f.write('TIME = {} hours since 1900-01-01 00:00:00 +00:00\n'.format(dt))

    #np.savetxt(pfid,np.flipud(p/0.01),fmt='%.3f')
    #np.savetxt(ufid,np.flipud(u),fmt='%.3f')
    #np.savetxt(vfid,np.flipud(v),fmt='%.3f')

         # close files
        for f in fi:
           f.close()
     
            
    def run(self,**kwargs):
        
        calc_dir = kwargs.get('rpath', './')
        bin_path = kwargs.get('solver', './')
        
        ncores = kwargs.get('ncores', './')
        
        
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
        
        
        # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./run_flow2d3d.sh {} {}'.format(bin_path,ncores)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
#        for line in iter(ex.stderr.readline,b''): print line
#        ex.stderr.close() 
        with open(calc_dir+'run.log', 'w') as f: 
          for line in iter(ex.stdout.readline,b''): 
            f.write(line)   
            sys.stdout.write(line)
            sys.stdout.flush()  
        ex.stdout.close()            
                 
#        exitCode = ex.returncode
#        print exitCode
#        if (exitCode == 0):
#        else:     
#            for line in iter(ex.stderr.readline,b''): 
#                sys.stdout.write(line)
#                sys.stdout.flush()
#            ex.stderr.close()        
                 
    
    def output(self,**kwargs):
        
         path = kwargs.get('rpath', './')
        
         with open(path+self.tag+'.pkl', 'w') as f:
               pickle.dump(self.__dict__,f)
        
    def save(self,**kwargs):
        
         path = kwargs.get('rpath', './')
        
         lista = [key for key, value in self.__dict__.items() if key not in ['meteo','dem']]
         dic = {k: self.__dict__.get(k, None) for k in lista}
         dic['model'] = self.__class__.__name__
         dic = dict(dic, **{k: self.__dict__.get(k, None).impl.__class__.__name__ for k in ['meteo','dem']})
         with open(path+'info.pkl', 'w') as f:
               pickle.dump(dic,f)
    
    @staticmethod
    def read(filename,**kwargs):
        
        filename = kwargs.get('filename', './')               
        with open(filename, 'r') as f:
              info=pickle.load(f)
        m = model(**info)
        return m
        
    def set(self,**kwargs):      
          
        path = kwargs.get('rpath', './') 
        
        #save mdf 
        mdf.write(self.mdf.inp, path+self.tag+'.mdf',selection=self.mdf.order)
        
        #save grid
        self.grid.write(path+self.tag+'.grd')
        
        #save dem
        # Write bathymetry file
        ba = Dep()
        # append the line/column of nodata 
        nodata=np.empty(self.ni)
        nodata.fill(np.nan)
        bat1=np.vstack((self.dem.impl.ival,nodata))
        nodata=np.empty((self.nj+1,1))
        nodata.fill(np.nan)
        bat2=np.hstack((bat1,nodata))
        ba.val = bat2
        ba.shape = bat2.shape

        Dep.write(ba,path+self.tag+'.dep')
        
        #save meteo
        self.force(rpath=path)
        
        #save bca
        
        
        #save obs
        
        
        #save enc
        #write enc out
        with open(path+self.tag+'.enc','w') as f:
            f.write('{:>5}{:>5}\n'.format(self.ni+1,1))  # add one like ddb
            f.write('{:>5}{:>5}\n'.format(self.ni+1,self.nj+1))
            f.write('{:>5}{:>5}\n'.format(1,self.nj+1))
            f.write('{:>5}{:>5}\n'.format(1,1))
            f.write('{:>5}{:>5}\n'.format(self.ni+1,1))
        
            
        
        
        