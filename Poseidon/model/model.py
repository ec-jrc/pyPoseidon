import os
import datetime
import numpy as np
import pickle
import xml.dom.minidom as md
from shutil import copy2
import subprocess
import sys
import pkg_resources

class model:
    def __init__(self, **kwargs):
        self.properties = kwargs.get('properties', {})
        
    def set(self):
        raise NotImplementedError("Subclass must implement abstract method")        
        
        
    def force(self):
        raise NotImplementedError("Subclass must implement abstract method")    
            

    def run(self):
        raise NotImplementedError("Subclass must implement abstract method")  
        

    def save(self):
        raise NotImplementedError("Subclass must implement abstract method")    
          
        
class d3d(model):
      
    def set(self,**kwargs):
        
        self.lon0 = kwargs.get('lon0', None)
        self.lon1 = kwargs.get('lon1', None)
        self.lat0 = kwargs.get('lat0', None)
        self.lat1 = kwargs.get('lat1', None)       
        self.date = kwargs.get('date', None)
        self.tag = kwargs.get('tag', None)
        self.ft1 = kwargs.get('ft1', None)
        self.ft2 = kwargs.get('ft2', None)     
        self.resolution = kwargs.get('resolution', None)
        self.x = kwargs.get('x', None)
        self.y = kwargs.get('y', None)     
        
                
        resmin=self.resolution*60
      
        # computei ni,nj / correct lat/lon

        if 'lon' in self.properties.keys():
            
          self.x=self.properties.lon
          self.y=self.properties.lat
              
          nj,ni=lon.shape
          self.lon0=lon.min()
          self.lon1=lon.max()
          self.lat0=lat.min()
          self.lat1=lat.max()

        else:

          ni=int(round((self.lon1-self.lon0)/self.resolution)) #these are cell numbers
          nj=int(round((self.lat1-self.lat0)/self.resolution))
  
          self.lon1=self.lon0+ni*self.resolution #adjust max lon to much the grid
          self.lat1=self.lat0+nj*self.resolution

       
        ni=ni+1 # transfrom to grid points
        nj=nj+1
        
        self.ni=ni
        self.nj=nj
        
        if self.x is None :
          sys.stdout.write('create grid')
          sys.stdout.flush()
          sys.stdout.write('\n')
        # set the grid 
          x=np.linspace(self.lon0,self.lon1,self.ni)
          y=np.linspace(self.lat0,self.lat1,self.nj)
          self.x,self.y=np.meshgrid(x,y)
        
    
    def force(self,**kwargs):
    
        path = kwargs.get('path', './')
        curvi = kwargs.get('curvi', False)
        
        dlat=self.lats[1,0]-self.lats[0,0]
        dlon=self.lons[0,1]-self.lons[0,0]
        lat0=self.lats[0,0] 
        lon0=self.lons[0,0] 

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
              f.write('n_cols           = {}\n'.format(self.u.shape[2]))
              f.write('n_rows           = {}\n'.format(self.u.shape[1]))
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
        for it in range(self.ft1,self.ft2): # nt + 0 hour
          ntime=self.date+datetime.timedelta(hours=it)
          dt=(ntime-time0).total_seconds()/3600.


          for f in fi:
             f.write('TIME = {} hours since 2000-01-01 00:00:00 +00:00\n'.format(dt))

          np.savetxt(pfid,np.flipud(self.p[it,:,:]/0.01),fmt='%.3f')
          np.savetxt(ufid,np.flipud(self.u[it,:,:]),fmt='%.3f')
          np.savetxt(vfid,np.flipud(self.v[it,:,:]),fmt='%.3f')

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
        
        calc_dir = kwargs.get('run_path', './')
        bin_path = kwargs.get('d3d', './')
        arch = kwargs.get('arch', './')
        
        ncores = kwargs.get('ncores', './')
        
        #retrieve the module path
        DATA_PATH = pkg_resources.resource_filename('Poseidon', 'misc/')
        
      # edit and save config file
        copy2(DATA_PATH + 'config_d_hydro.xml',calc_dir+'config_d_hydro.xml')

        xml=md.parse(calc_dir+'config_d_hydro.xml')

        xml.getElementsByTagName('mdfFile')[0].firstChild.replaceWholeText(self.tag+'.mdf')

        with open(calc_dir+'config_d_hydro.xml','w') as f:
            xml.writexml(f)

        copy2(DATA_PATH + 'run_flow2d3d.sh',calc_dir+'run_flow2d3d.sh')
        
        #make the script executable
        execf = calc_dir+'run_flow2d3d.sh'
        mode = os.stat(execf).st_mode
        mode |= (mode & 0o444) >> 2    # copy R bits to X
        os.chmod(execf, mode)
        
        
        # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./run_flow2d3d.sh {} {} {}'.format(bin_path,arch,ncores)], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
        for line in iter(ex.stderr.readline,b''): print line
        ex.stderr.close()
        for line in iter(ex.stdout.readline,b''): print line
        ex.stdout.close()    
    
    def save(self,**kwargs):
        
         path = kwargs.get('path', './')
        
         with open(path+self.tag+'.pkl', 'w') as f:
               pickle.dump(self.__dict__,f)
        
        
            
        