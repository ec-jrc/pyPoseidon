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

#local modules
import mdf
from grid import *
from dep import *
from bnd import *
from pyPoseidon.meteo import meteo
from pyPoseidon.dem import dem
from pyPoseidon.utils.get_value import get_value

#retrieve the module path
DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc/')
#DATA_PATH = os.path.dirname(Poseidon.__file__)+'/misc/'    
#info_data = ('lon0','lon1','lat0','lat1','date','tag','resolution','ft1','ft2')

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
        
        self.lon0 = kwargs.get('lon0', None)
        self.lon1 = kwargs.get('lon1', None)
        self.lat0 = kwargs.get('lat0', None)
        self.lat1 = kwargs.get('lat1', None)       
        self.date = kwargs.get('date', None)
        self.tag = kwargs.get('tag', None)
        self.resolution = kwargs.get('resolution', None)
        self.ft1 = kwargs.get('ft1', None)
        self.ft2 = kwargs.get('ft2', None)
        self.tide = kwargs.get('tide', False)
        self.atm = kwargs.get('atm', False)

        for attr, value in kwargs.iteritems():
                setattr(self, attr, value)
            
        self.model = self.__dict__.copy()    
        self.model['solver'] = self.__class__.__name__    
                      
    def set(self,**kwargs):

        gx = get_value(self,kwargs,'x',None)#kwargs.get('x', None)
        gy = get_value(self,kwargs,'y',None)#kwargs.get('y', None)    
        mdf_file = kwargs.get('mdf', None)  
        Tstart = self.date.hour*60+self.ft1*60     
        Tstop = self.date.hour*60+self.ft2*60
        step = get_value(self,kwargs,'step',None)#kwargs.get('step', None)
        rstep = get_value(self,kwargs,'rstep',None)#kwargs.get('rstep', None)
                                        
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

            if not self.atm: inp['Sub1'] = ' '

          # set tide only run
            if self.tide :
                inp['Filbnd']=self.tag+'.bnd'
                inp['Filana']=self.tag+'.bca'
 #           if 'Tidfor' not in order: order.append('Tidfor')
 #           inp['Tidfor']=[['M2','S2','N2','K2'], \
 #                       ['K1','O1','P1','Q1'], \
 #                         ['-----------']]
 
          # specify ini file
        # if 'Filic' not in order: order.append('Filic')
        # inp['Filic']=basename+'.ini'

          # netCDF output
        if 'FlNcdf' not in self.mdf.order: self.mdf.order.append('FlNcdf')

        self.mdf.inp['FlNcdf'] = 'map his'

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
        self.meteo = meteo(**z)  
        
        #dem
    def bath(self,**kwargs):
        z = self.__dict__.copy()        
        
        z['grid_x'] = self.grid.x
        z['grid_y'] = self.grid.y
        
        z.update(kwargs) 
        self.dem = dem(**z)  
        
    def bc(self,**kwargs):
        #define boundaries
        z = self.__dict__.copy()        
        
        z['lons'] = self.grid.x[0,:]
        z['lats'] = self.grid.y[:,0]
        
        ba = -self.dem.impl.ival
        ba[ba<0]=np.nan
        z['dem']=ba
        z['cn']=10
        
        z.update(kwargs) 
                
        self.bound = box(**z)

    def tidebc(self,**kwargs):
    
        self.tide = tide()
        for key,val in self.bound.__dict__.iteritems():
        
        # compute tide constituents
            tval = []
            if len(val) > 0. :                   
                blons=[]
                blats=[]
                for l1,l2 in val:
                    blons.append(self.grid.x[l1[1]-1,l1[0]-1])   
                    blats.append(self.grid.y[l1[1]-1,l1[0]-1])
                    blons.append(self.grid.x[l2[1]-1,l2[0]-1])   
                    blats.append(self.grid.y[l2[1]-1,l2[0]-1])
                       
                blons = np.array(blons)#.ravel().reshape(-1,2)[:,0]
                blats =  np.array(blats)#.ravel().reshape(-1,2)[:,1] 
            #                  print bound,blons,blats
                             
                tval = tide(tmodel=self.tmodel, tpath=self.tpath, blons=blons,blats=blats)
                    
            setattr(self.tide, key, tval)        
                                               
    
    def uvp(self,**kwargs):
                
        path = get_value(self,kwargs,'rpath','./') 
                    
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

          np.savetxt(pfid,np.flipud(self.meteo.impl.p[it,:,:]),fmt='%.3f')
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
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
        
        bin_path = get_value(self,kwargs,'exec', None)   
            
        ncores = get_value(self,kwargs,'ncores',1)
                
        if not os.path.exists( calc_dir+'config_d_hydro.xml') :
            
          print calc_dir, DATA_PATH
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
         z['date']=z['date'].isoformat()
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
        slevel = get_value(self,kwargs,'slevel',None) 
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        #save mdf 
        mdf.write(self.mdf.inp, path+self.tag+'.mdf',selection=self.mdf.order)

        #save grid
        self.grid.write(path+self.tag+'.grd')
        
        #save dem
        bat = -self.dem.impl.ival #reverse for the hydro run
      # bat[bat<-slevel]=np.nan #mask dry points
        
        # Write bathymetry file
        ba = Dep()
        # append the line/column of nodata 
        nodata=np.empty(self.ni)
        nodata.fill(np.nan)
        bat1=np.vstack((bat,nodata))
        nodata=np.empty((self.nj+1,1))
        nodata.fill(np.nan)
        bat2=np.hstack((bat1,nodata))
        ba.val = bat2
        ba.shape = bat2.shape

        Dep.write(ba,path+self.tag+'.dep')
        
        #save meteo
        if self.atm:
            self.uvp(**kwargs)
             
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
        
        
        #save enc
        #write enc out
        with open(path+self.tag+'.enc','w') as f:
            f.write('{:>5}{:>5}\n'.format(self.ni+1,1))  # add one like ddb
            f.write('{:>5}{:>5}\n'.format(self.ni+1,self.nj+1))
            f.write('{:>5}{:>5}\n'.format(1,self.nj+1))
            f.write('{:>5}{:>5}\n'.format(1,1))
            f.write('{:>5}{:>5}\n'.format(self.ni+1,1))
        
        
        
