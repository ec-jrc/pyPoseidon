"""
Schism model of pyPoseidon. It controls the creation, output & execution of a complete simulation based on schism
               
"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import os
import datetime
import numpy as np
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
import logging

#local modules
import pyPoseidon
import pyPoseidon.grid as pgrid
import pyPoseidon.meteo as pmeteo
import pyPoseidon.dem as pdem
from pyPoseidon.utils.get_value import get_value
from pyPoseidon.utils.converter import myconverter

#logging setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(levelname)-8s %(asctime)s:%(name)s:%(message)s')

file_handler = logging.FileHandler('model.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)

sformatter = logging.Formatter('%(levelname)-8s %(message)s')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(sformatter)

logger.addHandler(file_handler)
logger.addHandler(stream_handler)

#retrieve the module path
#DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc')
DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    
        
class schism():
     
    def __init__(self,**kwargs):
         
        rfolder = kwargs.get('rfolder', None)
        if rfolder:
            self.read_folder(**kwargs)
                
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
        
        if not hasattr(self, 'end_date'):
            #--------------------------------------------------------------------- 
            logger.warning('model not set properly, No end_date\n')
            #--------------------------------------------------------------------- 
            
        
        self.tag = kwargs.get('tag', 'schism')
        self.tide = kwargs.get('tide', False)
        self.atm = kwargs.get('atm', True)
    
        self.epath = kwargs.get('epath', None)
    
        self.solver = self.__class__.__name__    
        
                                                   
        for attr, value in kwargs.items():
            if not hasattr(self, attr): setattr(self, attr, value)  
        
        
#============================================================================================        
# CONFIG
#============================================================================================

    def config(self, config_file=None, output=False, **kwargs): 
        
        dic = get_value(self,kwargs,'parameters',None)
#        param_file = get_value(self,kwargs,'config_file',None)
           
        if config_file :
            #--------------------------------------------------------------------- 
            logger.info('reading parameter file {}\n'.format(config_file))
            #---------------------------------------------------------------------             
        else:
            #--------------------------------------------------------------------- 
            logger.info('using default parameter file ...\n')
            #---------------------------------------------------------------------             
            
            config_file = DATA_PATH + 'param.in.sample'
        
        params = pd.read_csv(config_file,engine='python',comment='!', header=None, delimiter=' =', names=['attrs','vals'])

        fix_list = [k for k in params.attrs if '=' in k ]

        for k in fix_list:
            try:
                name, value = params.loc[params.attrs == k,'attrs'].values[0].split('=')
                params.loc[params.attrs == k,'attrs'] = name
                params.loc[params.attrs == name,'vals'] = value
            except:
                pass
            
        params = params.set_index('attrs')
        params.index = params.index.str.strip()
        params.vals = params.vals.str.strip()
        
        #update key values
        
        params.loc['start_year'] = self.start_date.year
        params.loc['start_month'] = self.start_date.month
        params.loc['start_day'] = self.start_date.day
        params.loc['start_hour'] = self.start_date.hour
                
        #update 
        if dic :
            for key,val in dic.items():
                params.loc[key] = val
        
        #test rnday
        if np.float(params.loc['rnday'].values[0]) * 24 * 3600 > (self.end_date - self.start_date).total_seconds():
            #--------------------------------------------------------------------- 
            logger.warning('rnday larger than simulation range\n')
            logger.warning('rnday={} while simulation time is {}\n'.format(params.loc['rnday'].values[0],(self.end_date - self.start_date).total_seconds()/(3600*24.)))
            #---------------------------------------------------------------------             
            
        
        self.params = params
        
                        
        if output: 
            #save params 
            #--------------------------------------------------------------------- 
            logger.info('output param.in file ...\n')
            #---------------------------------------------------------------------             
            
            path = get_value(self,kwargs,'rpath','./') 
            self.params.to_csv(path + 'param.in', header=None, sep='=')  #save to file

#============================================================================================        
# METEO
#============================================================================================       
    def force(self,**kwargs):
        
        meteo_files =  get_value(self,kwargs,'meteo_files',None)        
        
        kwargs.update({'meteo_files':meteo_files})
        
        flag = get_value(self,kwargs,'update',[])
        
        z = {**self.__dict__, **kwargs} # merge self and possible kwargs
        
        # check if files exist
        if flag :     
            if ('meteo' in flag) | ('all' in flag) : 
                self.meteo = pmeteo.meteo(**z)        
            else:
                logger.info('skipping meteo ..\n')
        else:
            self.meteo = pmeteo.meteo(**z)
        
        # add 1 hour for Schism issue with end time   
        ap = self.meteo.Dataset.isel(time = -1)
        ap['time'] = ap.time.values + pd.to_timedelta('1H')
        
        self.meteo.Dataset = xr.concat([self.meteo.Dataset,ap],dim='time')
        

    @staticmethod 
    def to_force(ar0,**kwargs):
        
        logger.info('writing meteo files ..\n')
                
        path = kwargs.get('rpath','./') 
                
        [p,u,v] = kwargs.get('vars','[None,None,None]')                
        
        ar = ar0.sortby('latitude', ascending=True)
            
        xx, yy = np.meshgrid(ar.longitude.data, ar.latitude.data) 
        
        zero = np.zeros(ar[p].data.shape)
                       
        date = kwargs.get('date',ar.time[0].data)
        
        udate = pd.to_datetime(date).strftime('%Y-%m-%d')
                       
        bdate = pd.to_datetime(date).strftime('%Y %m %d %H').split(' ')
        
        tlist = (ar.time.data - pd.to_datetime([udate]).values).astype('timedelta64[s]')/3600.
        
        tlist = tlist.astype(float)/24.
        
        
        bdate = [int(q) for q in bdate[:3]] + [0]
        
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
        
        filename = kwargs.get('filename','sflux/sflux_air_1.001.nc') 
               
        sout.to_netcdf(path+filename)
                
#============================================================================================        
# DEM
#============================================================================================       

    def bath(self,**kwargs):
 #       z = self.__dict__.copy()        
        
        kwargs['grid_x'] = self.grid.Dataset.SCHISM_hgrid_node_x.values
        kwargs['grid_y'] = self.grid.Dataset.SCHISM_hgrid_node_y.values
        
        dpath =  get_value(self,kwargs,'dem_file',None)        
        
        kwargs.update({'dem_file':dpath})
                
        flag = get_value(self,kwargs,'update',[])
        # check if files exist
        if flag :
            if ('dem' in flag) | ('all' in flag) :
                self.dem = pdem.dem(**kwargs)
            else:
                logger.info('dem from grid file\n')
#============================================================================================        
# EXECUTION
#============================================================================================
    def create(self,**kwargs):

        if not kwargs : kwargs = self.__dict__.copy()
                                         
        # Grid         
        self.grid=pgrid.grid(type='tri2d',**kwargs)
                 
        # set lat/lon from file
        if hasattr(self, 'grid_file'):
            kwargs.update({'minlon' : self.grid.Dataset.SCHISM_hgrid_node_x.values.min()})
            kwargs.update({'maxlon' : self.grid.Dataset.SCHISM_hgrid_node_x.values.max()})
            kwargs.update({'minlat' : self.grid.Dataset.SCHISM_hgrid_node_y.values.min()})
            kwargs.update({'maxlat' : self.grid.Dataset.SCHISM_hgrid_node_y.values.max()})
                                     
        # get bathymetry
        self.bath(**kwargs)

        # get boundaries
        # self.bc()
                
        #get meteo
        if self.atm :  self.force(**kwargs)
        
        #get tide
        if self.tide : self.tidebc()
        
        
        self.config(**kwargs)


    def output(self,**kwargs):      
        
        path = get_value(self,kwargs,'rpath','./') 
        flag = get_value(self,kwargs,'update',[])
        
        
        if not os.path.exists(path):
            os.makedirs(path)
        
        # save sflux_inputs.txt
        if not os.path.exists(path+'sflux'):
            os.makedirs(path+'sflux')
        
        with open(path+'sflux/sflux_inputs.txt', 'w') as f:
            f.write('&sflux_inputs\n')
            f.write('/ \n\n')
            
        # save bctides.in
        nobs = [key for key in self.grid.Dataset.keys() if 'open' in key]
        
        with open(path + 'bctides.in', 'w') as f:
            f.write('Header\n')
            f.write('{} {}\n'.format(0, 40.)) #  ntip tip_dp
            f.write('{}\n'.format(0)) #nbfr
            f.write('{}\n'.format(len(nobs))) #number of open boundaries
            for i in range(len(nobs)):
                nnodes = self.grid.Dataset[nobs[i]].dropna(dim='dim_0').size
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
        
            bat = -self.dem.Dataset.ival.values.astype(float) #minus for the hydro run
                            
            self.grid.Dataset.depth.loc[:bat.size] = bat
                            
            self.grid.to_file(filename= path+'hgrid.gr3')
            copyfile(path+'hgrid.gr3', path+'hgrid.ll')
    
            logger.info('updating bathymetry ..\n')
        
        except AttributeError as e:
            
            logger.info('Keeping bathymetry from hgrid.gr3 ..\n')
        
            copyfile(self.grid_file, path+'hgrid.gr3') #copy original grid file
            copyfile(path+'hgrid.gr3', path+'hgrid.ll')                
                 
                 
        # manning file
        manfile=path+'manning.gr3'
        
        
        if hasattr(self, 'manfile') :
            copyfile(self.manfile, manfile) #copy original grid file
            if self.manfile == manfile:
                logger.info('Keeping manning file ..\n')
        
        
        if hasattr(self, 'manning') :
            nn = self.grid.Dataset.nSCHISM_hgrid_node.size
            n3e = self.grid.Dataset.nSCHISM_hgrid_face.size
        
            with open(manfile,'w') as f:
                f.write('\t 0 \n')
                f.write('\t {} {}\n'.format(n3e,nn))
                
            df = self.grid.Dataset[['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']].to_dataframe()
                
            df['man'] = self.manning
            
            df.index = np.arange(1, len(df) + 1)
                
            df.to_csv(manfile,index=True, sep='\t', header=None,mode='a', float_format='%.10f',columns=['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','man'] )
            
            logger.info('Manning file created..\n')
                
        
        # windrot_geo2proj
        
        windfile=path+'windrot_geo2proj.gr3'

        if hasattr(self, 'windproj') :
            copyfile(self.windproj, windfile) #copy original grid file
            if self.windproj != windproj :
                logger.info('Keeping windproj file ..\n')

        if hasattr(self, 'windrot') :
            with open(windfile,'w') as f:
                f.write('\t 0 \n')
                f.write('\t {} {}\n'.format(n3e,nn))
                    
            df = self.grid.Dataset[['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','depth']].to_dataframe()
               
            df['windrot'] = self.windrot
            
            df.index = np.arange(1, len(df) + 1)
        
            df.to_csv(windfile,index=True, sep='\t', header=None,mode='a', float_format='%.10f',columns=['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y','windrot'] )
                
            logger.info('Windrot_geo2proj file created..\n')
            
        #save meteo
        if hasattr(self, 'atm') :
           try:
              self.to_force(self.meteo.Dataset,vars=['msl','u10','v10'],rpath=path,**kwargs)
           except AttributeError as e:
              logger.warning('no meteo data available.. no update..\n')
              pass
                              
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
                        
        bin_path = get_value(self,kwargs,'epath', None) 
        
        if bin_path is None:
            #------------------------------------------------------------------------------ 
            logger.warning('Schism executable path (epath) not given -> using default \n')
            #------------------------------------------------------------------------------
            bin_path = 'schism'
              
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

        #--------------------------------------------------------------------- 
        logger.info('output done\n')
        #--------------------------------------------------------------------- 

                                         
                
    def run(self,**kwargs):
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
                
        bin_path = get_value(self,kwargs,'epath', None)   
            
        ncores = get_value(self,kwargs,'ncores',1)
        
        conda_env = get_value(self,kwargs,'conda_env', None)

        #--------------------------------------------------------------------- 
        logger.info('executing model\n')
        #--------------------------------------------------------------------- 

            # note that cwd is the folder where the executable is
        ex=subprocess.Popen(args=['./launchSchism.sh'], cwd=calc_dir, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=1)
            
        with open(calc_dir+'err.log', 'w') as f: 
          for line in iter(ex.stderr.readline,b''): 
            f.write(line.decode(sys.stdout.encoding))   
            logger.info(line.decode(sys.stdout.encoding))
        ex.stderr.close()            

        with open(calc_dir+'run.log', 'w') as f: 
          for line in iter(ex.stdout.readline,b''): 
            f.write(line.decode(sys.stdout.encoding))   
            logger.info(line.decode(sys.stdout.encoding))
        ex.stdout.close()         
        
        #--------------------------------------------------------------------- 
        logger.info('FINISHED\n')
        #--------------------------------------------------------------------- 
        
                                  
    def save(self,**kwargs):
               
         path = get_value(self,kwargs,'rpath','./')
        
         lista = [key for key, value in self.__dict__.items() if key not in ['meteo','dem','grid']]
         dic = {k: self.__dict__.get(k, None) for k in lista}

         grid=self.__dict__.get('grid', None)
         if isinstance(grid,np.str):
             dic.update({'grid':grid})
         else:
             dic.update({'grid':grid.__class__.__name__})
         
         dem=self.__dict__.get('dem', None)
         if isinstance(dem,np.str):
             dic.update({'dem':dem})
         elif isinstance(dem,pdem.dem):
             dic.update({'dem': dem.Dataset.elevation.attrs})

         meteo=self.__dict__.get('meteo', None)
         if isinstance(meteo,np.str):
             dic.update({'meteo':meteo})
         elif isinstance(meteo,pmeteo.meteo):
             dic.update({'meteo':meteo.Dataset.attrs})

         dic['version']=pyPoseidon.__version__
               
         for attr, value in dic.items():
             if isinstance(value, datetime.datetime) : dic[attr]=dic[attr].isoformat()
             if isinstance(value, pd.Timedelta) : dic[attr]=dic[attr].isoformat()          
             if isinstance(value, pd.DataFrame) : dic[attr]=dic[attr].to_dict()
         json.dump(dic,open(path+self.tag+'_model.json','w'),default = myconverter)     
         
         
    def execute(self,**kwargs):
        
        self.create(**kwargs)
        self.output(**kwargs) 
        self.save(**kwargs)
        self.run(**kwargs)
                 
    def read_folder(self, rfolder,**kwargs):
        
        s = glob.glob(rfolder + '/param*')
        mfiles = glob.glob(rfolder + '/sflux/*.nc')
        
        hfile = rfolder + '/hgrid.gr3' # Grid
        params = pd.read_csv(s[0],engine='python',comment='!', header=None, delimiter=' =', names=['attrs','vals'])

        fix_list = [k for k in params.attrs if '=' in k ]

        for k in fix_list:
            try:
                name, value = params.loc[params.attrs == k,'attrs'].values[0].split('=')
                params.loc[params.attrs == k,'attrs'] = name
                params.loc[params.attrs == name,'vals'] = value
            except:
                pass
            
        self.params = params.set_index('attrs')
        self.grid = pgrid.grid(type='tri2d',grid_file=hfile)
        
        self.meteo = xr.open_mfdataset(mfiles) # Meteo
        
