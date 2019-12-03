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
import geopandas as gp
import glob
from shutil import copyfile
import xarray as xr
import cartopy.feature as cf
import geopandas as gp


#local modules
import pyPoseidon
import pyPoseidon.grid as pgrid
import pyPoseidon.meteo as pmeteo
import pyPoseidon.dem as pdem
from pyPoseidon.utils.get_value import get_value
from pyPoseidon.utils.converter import myconverter
from pyPoseidon.utils import obs
from pyPoseidon.utils.cpoint import closest_node

import logging
logger = logging.getLogger('pyPoseidon')

#retrieve the module path
#DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc')
DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    
        
class schism():
     
    def __init__(self,**kwargs):
         
        rfolder = kwargs.get('rfolder', None)
        if rfolder:
            self.read_folder(**kwargs)

        self.geometry = kwargs.get('geometry', None)
        
        if self.geometry:
                
            if isinstance(self.geometry,dict): 
                self.lon_min = self.geometry['lon_min']
                self.lon_max = self.geometry['lon_max']
                self.lat_min = self.geometry['lat_min']
                self.lat_max = self.geometry['lat_max']
            elif isinstance(self.geometry,str):
                  
                try:
                    geo = gp.GeoDataFrame.from_file(self.geometry)
                except:
                    logger.error('geometry argument not a valid geopandas file')
                    sys.exit(1)
            
                self.lon_min, self.lat_min, self.lon_max, self.lat_max = geo.total_bounds


        # coastlines
        coastlines = kwargs.get('coastlines',None)
        
        if coastlines is None:
            cr = kwargs.get('coast_resolution', 'l')
    
            # world polygons - user input
            coast = cf.NaturalEarthFeature(
                category='physical',
                name='land',
                scale='{}m'.format({'l':110, 'i':50, 'h':10}[cr]))
    
           
            self.coastlines = gp.GeoDataFrame(geometry = [x for x in coast.geometries()])    
    
        else:
            
            try:
                coast = gp.GeoDataFrame.from_file(coastlines)
            except:
                coast = gp.GeoDataFrame(coastlines)
                
            self.coastlines = coast
        
               
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
        self.monitor = kwargs.get('monitor', False)
    
        try:
            self.epath = os.environ['SCHISM']
        except:
            self.epath = kwargs.get('epath', None)
    
        self.solver = self.__class__.__name__    
                                                       
        for attr, value in kwargs.items():
            if not hasattr(self, attr): setattr(self, attr, value)  
        
        setattr(self, 'misc', {})
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
                
        meteo_source =  get_value(self,kwargs,'meteo_source',None)        
        
        kwargs.update({'meteo_source':meteo_source})
        
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
        
        if hasattr(self, 'meteo'):
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
        
        dpath =  get_value(self,kwargs,'dem_source',None)        
        
        kwargs.update({'dem_source':dpath})
                
        flag = get_value(self,kwargs,'update',[])
        # check if files exist
        if flag :
            if ('dem' in flag) | ('all' in flag) :
                kwargs.update({'lon_min': self.lon_min, 'lat_min': self.lat_min, 'lon_max':self.lon_max, 'lat_max': self.lat_max})
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
            kwargs.update({'lon_min' : self.grid.Dataset.SCHISM_hgrid_node_x.values.min()})
            kwargs.update({'lon_max' : self.grid.Dataset.SCHISM_hgrid_node_x.values.max()})
            kwargs.update({'lat_min' : self.grid.Dataset.SCHISM_hgrid_node_y.values.min()})
            kwargs.update({'lat_max' : self.grid.Dataset.SCHISM_hgrid_node_y.values.max()})
            
            self.lon_min = self.grid.Dataset.SCHISM_hgrid_node_x.values.min()
            self.lon_max = self.grid.Dataset.SCHISM_hgrid_node_x.values.max()
            self.lat_min = self.grid.Dataset.SCHISM_hgrid_node_y.values.min()
            self.lat_max = self.grid.Dataset.SCHISM_hgrid_node_y.values.max()
            
                                     
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
        split_by = get_value(self,kwargs,'split_by',None)
        
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
                nnodes = self.grid.Dataset[nobs[i]].dropna(dim='index').size
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
              if split_by :
                  times, datasets = zip(*self.meteo.Dataset.groupby('time.{}'.format(split_by)))
                  mpaths = ['sflux/sflux_air_1.{:03d}.nc'.format(t + 1) for t in np.arange(len(times))]
                  for das,mpath in list(zip(datasets,mpaths)):
                      self.to_force(das,vars=['msl','u10','v10'],rpath=path, filename=mpath, **kwargs)
              else:  
                  self.to_force(self.meteo.Dataset,vars=['msl','u10','v10'],rpath=path,**kwargs)                 
           except AttributeError as e:
              logger.warning('no meteo data available.. no update..\n')
              pass
                              
        
        calc_dir = get_value(self,kwargs,'rpath','./') 
        
        try:
            bin_path = os.environ['SCHISM']
        except:
            bin_path = get_value(self,kwargs,'epath', None)
                                
        if bin_path is None:
            #------------------------------------------------------------------------------ 
            logger.warning('Schism executable path (epath) not given -> using default \n')
            #------------------------------------------------------------------------------
            bin_path = 'schism'
              
        ncores = get_value(self,kwargs,'ncores',1)
                            
            
        with open(calc_dir + 'launchSchism.sh', 'w') as f:
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
                                    
        ncores = get_value(self,kwargs,'ncores',1)
        
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
        if self.monitor : self.set_obs() 
        self.save(**kwargs)
        self.run(**kwargs)
                 
    def read_folder(self, rfolder,**kwargs):
        
        self.rpath = rfolder
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
        
        #meteo 
        ma = []
        for ifile in mfiles:
            g = xr.open_dataset(ifile)
            ts = '-'.join(g.time.attrs['base_date'].astype(str)[:3])
            time_r = pd.to_datetime(ts) 
            times = time_r + pd.to_timedelta(g.time.values,unit='D').round('H')
            g = g.assign_coords({'time':times})
            ma.append(g)
        
        try:
            self.meteo = xr.concat(ma,dim='time')
        except:
            logger.warning('No meteo files loaded')
            pass

    def global2local(self,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
                    
        #Read the global node index distribution to the cores
        gfiles = glob.glob(path+'outputs/local_to_global_*')
        gfiles.sort()
        
        #create a dict from filenames to identify parts in the dataframes below
        keys=[]
        for name in gfiles:
            keys.append('core{}'.format(name.split('/')[-1].split('_')[-1]))

        # Parsing the files
        l2g=[]
        for i in range(len(gfiles)):
            l2g.append(pd.read_csv(gfiles[i],header=None,delim_whitespace=True,engine='python'))

        # We read from the first file the header (it is the same for all)
        header = l2g[0].iloc[0]
        header.index = ['ns_global','ne_global','np_global','nvrt','nproc','ntracers','T','S','GEN','AGE','SED3D','EcoSim','ICM','CoSINE','Feco','TIMOR','FABM']

        #get the number of elems from all files
        nels = []
        for i in range(len(l2g)):
            nels.append(int(l2g[i].iloc[2].dropna()[0]))

        frames = []
        for i in range(len(l2g)):
            df = l2g[i].iloc[3:nels[i]+3,l2g[i].columns[:2]].astype(int)
            df = df.reset_index(drop=True)
            frames.append(df)

        elems = pd.concat(frames,keys=keys)

        elems.columns = ['local','global_n']

        #get the number of nodes from all files
        nq= []
        for i in range(len(l2g)):
            nq.append(int(l2g[i].iloc[nels[i]+3].dropna()[0]))

        nframes = []
        for i in range(len(l2g)):
            df = l2g[i].iloc[nels[i]+4:nels[i]+4+nq[i],l2g[i].columns[:2]].astype(int)
            df = df.reset_index(drop=True)
            nframes.append(df)

        nodes = pd.concat(nframes,keys=keys)

        nodes.columns = ['local','global_n']

        #get the number of nodes from all files
        nw= []
        for i in range(len(l2g)):
            j = nels[i]+4+nq[i]
            nw.append(int(l2g[i].iloc[j].dropna()[0]))
            
        wframes = []
        for i in range(len(l2g)):
            j0 = nels[i]+5+nq[i]
            j1 = j0+nw[i]
            df = l2g[i].iloc[j0:j1,l2g[i].columns[:2]].astype(int)
            df = df.reset_index(drop=True)
            wframes.append(df)

        re = pd.concat(wframes,keys=keys)

        re.columns = ['local','global_n']

        # Read secondary headers
        h0 = l2g[0].iloc[nels[0] + nq[0] + nw[0] + 6].dropna()
        h0.index = ['start_year','start_month','start_day','start_hour','utc_start']

        h1 = l2g[0].iloc[nels[0] + nq[0] + nw[0] + 7].dropna()
        h1.index = ['nrec','dtout','nspool','nvrt','kz','h0','h_s','h_c','theta_b','theta_f','ics']
        h1 = h1.apply(pd.to_numeric)

        ztots = ['ztot_'+str(i) for i in range(1,h1.kz.astype(int)-1)]

        sigmas = ['sigma_'+str(i) for i in range(h1.nvrt.astype(int) - h1.kz.astype(int) + 1) ]

        h2 = l2g[0].iloc[nels[0] + nq[0] + nw[0] + 8].dropna()
        h2.index = ztots + sigmas

        #combine headers
        self.misc.update({'header' : pd.concat([h0, h1, h2])})

        #Read grid
        gframes = []
        for i in range(len(l2g)):
            j0 = nels[i] + nq[i] + nw[i] + 10
            j1 = j0+nq[i]
            df = l2g[i].iloc[j0:j1,l2g[i].columns[:4]].astype(float)
            df = df.reset_index(drop=True)
            gframes.append(df)

        grid = pd.concat(gframes,keys=keys)

        grid.columns = ['lon','lat','depth','kbp00']

        #Droping duplicates
        cnodes = nodes.global_n.drop_duplicates() # drop duplicate global nodes and store the values to an array
        grid = grid.drop_duplicates() # keep only one of the duplicates (the first by default)
        grid.index = grid.index.droplevel() # drop multi-index
        grid = grid.reset_index(drop=True) # reset index
        grid.index = cnodes.values - 1 # reindex based on the global index, -1 for the python convention
        grd = grid.sort_index() #sort with the new index (that is the global_n)
        self.misc.update({'grd' : grd.reset_index(drop=True)})#reindex for final version

        #Read tessalation
        eframes = []
        for i in range(len(l2g)):
            j0 = nels[i] + nq[i] + nw[i] + 10 + nq[i]
            j1 = j0+nels[i]
            df = l2g[i].iloc[j0:j1,l2g[i].columns[:4]].astype(int)
            df = df.reset_index(drop=True)
            eframes.append(df)

        tri = pd.concat(eframes,keys=keys)

        tri.columns = ['type','a','b','c']

        # Replace local index with the global one 
        for key in keys:
            nod = nodes.loc[key].copy() # make a copy of the the core
            nod = nod.set_index('local') # reset the index using the local values (in essense set the index to start from 1...)
            tri.loc[key,'ga'] = nod.reindex(tri.loc[key,'a'].values).values
            tri.loc[key,'gb'] = nod.reindex(tri.loc[key,'b'].values).values
            tri.loc[key,'gc'] = nod.reindex(tri.loc[key,'c'].values).values 
        
        tri.loc[:,'ga'] = tri.loc[:,'ga'].apply(pd.to_numeric(int)) # make integer
        tri.loc[:,'gb'] = tri.loc[:,'gb'].apply(pd.to_numeric(int)) # make integer
        tri.loc[:,'gc'] = tri.loc[:,'gc'].apply(pd.to_numeric(int)) # make integer
        
        tri = tri.drop(['a','b','c'],axis=1) #drop local references

        #sort
        gt3 = tri.loc[:,['ga','gb','gc']].copy() # make a copy
        gt3.index = gt3.index.droplevel() # drop multi-index
        gt3 = gt3.reset_index(drop=True)
        # Now we need to put them in order based on the global index in elems
        gt3.index = elems.global_n.values # we set the index equal to the global_n column
        gt3 = gt3.sort_index() #sort them
        
        #add nan column in place of the fourth node. NOTE:  This needs to be tested for quadrilaterals
        gt3['gd']=np.nan
        
        gt3 = gt3.reset_index() # reset to add more columns without problems
        
        ## Add mean x, y of the elememts. To be used in the output
        gt3['x1'] = grd.loc[gt3['ga'].values - 1, 'lon'].values #lon of the index, -1 for python convention
        gt3['y1'] = grd.loc[gt3['ga'].values - 1, 'lat'].values #lat of the index
        gt3['x2'] = grd.loc[gt3['gb'].values - 1, 'lon'].values
        gt3['y2'] = grd.loc[gt3['gb'].values - 1, 'lat'].values
        gt3['x3'] = grd.loc[gt3['gc'].values - 1, 'lon'].values
        gt3['y3'] = grd.loc[gt3['gc'].values - 1, 'lat'].values


        gt3['xc'] =  gt3[['x1', 'x2', 'x3']].mean(axis=1) #mean lon of the element
        gt3['yc'] =  gt3[['y1', 'y2', 'y3']].mean(axis=1)

        ## min kbe
        gt3['kbe1'] = grd.loc[gt3['ga'] - 1,'kbp00'].values
        gt3['kbe2'] = grd.loc[gt3['gb'] - 1,'kbp00'].values
        gt3['kbe3'] = grd.loc[gt3['gc'] - 1,'kbp00'].values
        #gt3['kbe4'] = grd.loc[gt3['gd'],'kbp00'].values

        gt3['kbe'] = gt3[['kbe1', 'kbe2', 'kbe3']].min(axis=1)
        
        self.misc.update({'gt3' : gt3.set_index('index')}) # set index back 
        
        
        
        #Droping duplicates
        self.misc.update({'melems' : elems.loc[elems.global_n.drop_duplicates().index]}) # create the retaining mask
        self.misc.update({'msides' : re.loc[re.global_n.drop_duplicates().index]}) # keep only one of the duplicates
        self.misc.update({'mnodes' : nodes.loc[nodes.global_n.drop_duplicates().index]}) # keep only one of the duplicates
                
    
    def hotstart(self, it=None, **kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
                          
        if not 'melems' in self.misc: self.global2local(**kwargs)

        hfiles = glob.glob(path+'outputs/hotstart_*_{}.nc'.format(it))
        hfiles.sort()

        #store them in a list
        out=[]
        for i in range(len(hfiles)):
            out.append(xr.open_dataset(hfiles[i]))

        # Create dataset      
        side=[]
        node=[]
        el=[]
        one=[]

        for key in out[0].variables:
            if 'nResident_side' in out[0][key].dims : 
                r = self.combine_(key,out,self.misc['msides'],'nResident_side')       
                side.append(r)
            elif 'nResident_node' in out[0][key].dims : 
                r = self.combine_(key,out,self.misc['mnodes'],'nResident_node')
                node.append(r)
            elif 'nResident_elem' in out[0][key].dims : 
                r = self.combine_(key,out,self.misc['melems'],'nResident_elem')
                el.append(r)
            elif len(out[0][key].dims) == 1:
                one.append(out[0][key])
        
        
        side = xr.merge(side).rename({'nResident_side':'side'})
        el = xr.merge(el).rename({'nResident_elem':'elem'})
        node = xr.merge(node).rename({'nResident_node':'node'})
        one = xr.merge(one).rename({'one':'one_new','it':'iths'})
        
        #merge
        xdat = xr.merge([side,el,node,one])
                            
        hfile = 'hotstart_it={}.nc'.format(xdat.iths.values[0])
        logger.info('saving hotstart file\n')

        xdat.to_netcdf(path + 'outputs/{}'.format(hfile))



    ## Any variable
    def combine(self, out, g2l, name):
        keys = g2l.index.get_level_values(0).unique()
        r=[]
        for i in range(len(out)):    
            v = out[i].to_pandas()
            v.index += 1   
            mask = g2l.loc[keys[i],'local']
            vm = v.loc[mask]
            vm.index = g2l.loc[keys[i],'global_n'].values
            r.append(vm)
        r = pd.concat(r).sort_index()
        r.index -= 1
        r.index.name = name
        return r
    
    def combine_(self, var, out, g2l, name):
        if len(out[0][var].shape) == 3:
            dd = []
            for i in range(out[0][var].shape[2]):
                o = []
                for j in range(len(out)):
                    o.append(out[j][var].loc[{out[j][var].dims[2]:i}])
                r = self.combine(o, g2l, name)
                dd.append(xr.DataArray(r.values, dims = out[0][var].dims[:2], name=var))
            
            tr = xr.concat(dd,dim=out[0][var].dims[2])
            a,b,c = out[0][var].dims
            tr = tr.transpose(a,b,c)
            return tr
    
        elif len(out[0][var].shape) < 3:
            o = []
            for j in range(len(out)):
                o.append(out[j][var])
            r = self.combine(o, g2l, name)
            return(xr.DataArray(r.values, dims = list(out[0][var].dims), name=var))



    def xcombine(self, tfs):

        # Create dataset      
        side=[]
        node=[]
        el=[]
        single=[]

        for key in tfs[0].variables:
            if 'nSCHISM_hgrid_face' in tfs[0][key].dims : 
                r = self.combine_(key,tfs,self.misc['melems'],'nSCHISM_hgrid_face')       
                side.append(r)
            elif 'nSCHISM_hgrid_node' in tfs[0][key].dims : 
                r = self.combine_(key,tfs,self.misc['mnodes'],'nSCHISM_hgrid_node')
                node.append(r)
            elif len(tfs[0][key].dims) == 1:
                single.append(tfs[0][key])


        side = xr.merge(side)
        el = xr.merge(el)
        node = xr.merge(node)
        single = xr.merge(single)

        #merge
        return xr.merge([side,el,node,single])
        
    
    def tcombine(self, hfiles, sdate, times):
    
        xall=[]
        for k in range(times.shape[0]):
            tfs=[]
            for i in range(len(hfiles)):
                tfs.append(xr.open_dataset(hfiles[i]).isel(time=k))

            xall.append(self.xcombine(tfs))
    
        xdat = xr.concat(xall,dim='time')
        xdat = xdat.assign_coords(time=times)

        return xdat
        
        
        
    
    #https://stackoverflow.com/questions/41164630/pythonic-way-of-removing-reversed-duplicates-in-list
    @staticmethod
    def remove_reversed_duplicates(iterable):
        # Create a set for already seen elements
        seen = set()
        for item in iterable:
            # Lists are mutable so we need tuples for the set-operations.
            tup = tuple(item)
            if tup not in seen:
                # If the tuple is not in the set append it in REVERSED order.
                seen.add(tup[::-1])
                # If you also want to remove normal duplicates uncomment the next line
                # seen.add(tup)
                yield item
    
    
    def read_vgrid(self,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./')
        
        vgrid = pd.read_csv(path + 'vgrid.in', header=None)
        
        self.misc.update({'ivcor' : vgrid.iloc[0].astype(int).values[0]})
    
        [Nz, kz, hs] = vgrid.iloc[1].str.split(' ')[0]
        
        self.misc.update({'Nz' : int(Nz)})
        self.misc.update({'kz' : int(kz)})
        self.misc.update({'hs' : float(hs)})
        
        zlevels = vgrid.iloc[3:3+self.misc['kz'],0].str.split(' ', n = 2, expand = True)
        zlevels.columns = ['level_index','z-coordinates']
        zlevels.set_index('level_index', inplace=True)

        self.misc.update({'zlevels' : zlevels})
        
        constants_index = 3+self.misc['kz']+1

        [h_c, theta_b, theta_f] = vgrid.iloc[constants_index].str.split(' ')[0]
        self.misc.update({'h_c' : float(h_c)})
        self.misc.update({'theta_b' : float(theta_b)})
        self.misc.update({'theta_f' : float(theta_f)})

        
        sl_index0 = constants_index+1
        sl_index1 = sl_index0 + self.misc['Nz'] - self.misc['kz'] + 1
        
        slevels = vgrid.iloc[sl_index0:sl_index1,0].str.split(' ', n = 2, expand = True)
        slevels.columns = ['level_index','s-coordinates']
        slevels.set_index('level_index', inplace=True)

        self.misc.update({'slevels' : slevels})
        
        
        
    def results(self,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
        
        if not hasattr(self, 'melems'): 
            logger.info('retrieving index references ... \n')
            self.global2local(**kwargs)
            logger.info('... done \n')
            
        # Create grid xarray Dataset
        grd = self.misc['grd']
        gt3 = self.misc['gt3']
        
        # node based variables
        grd.kbp00 = grd.kbp00.astype(int)
        xnodes = grd.to_xarray().rename({'lon':'SCHISM_hgrid_node_x','lat':'SCHISM_hgrid_node_y','kbp00':'node_bottom_index', 'index':'nSCHISM_hgrid_node'})

        xnodes = xnodes.drop_vars('nSCHISM_hgrid_node')

        # element based variables
        gt34 = gt3.loc[:,['ga','gb','gc']].values # SCHISM_hgrid_face_nodes
        xelems = xr.Dataset({
                         u'SCHISM_hgrid_face_nodes' : ([u'nSCHISM_hgrid_face', u'nMaxSCHISM_hgrid_face_nodes'], gt34 - 1),  #-> start index = 0 
                         u'SCHISM_hgrid_face_x' : ([u'nSCHISM_hgrid_face'], gt3.loc[:,'xc'].values),
                         u'SCHISM_hgrid_face_y' : ([u'nSCHISM_hgrid_face'], gt3.loc[:,'yc'].values),                 
                         u'ele_bottom_index': ([u'nSCHISM_hgrid_face'], gt3.kbe.values )})

        logger.info('done with node based variables \n')
        
        # edge based variables
        sides=[]
        for [ga,gb,gc] in gt3.loc[:,['ga','gb','gc']].values:
            sides.append([gb,gc])
            sides.append([gc,ga])
            sides.append([ga,gb])

        #removing duplicates
        sides = list(self.remove_reversed_duplicates(sides))    
            
        ed = pd.DataFrame(sides, columns=['node1','node2'])

        #mean x, y 
        ed['x1'] = grd.loc[ed['node1'].values - 1, 'lon'].values #lon of the index, -1 for python convention
        ed['y1'] = grd.loc[ed['node1'].values - 1, 'lat'].values #lat of the index
        ed['x2'] = grd.loc[ed['node2'].values - 1, 'lon'].values
        ed['y2'] = grd.loc[ed['node2'].values - 1, 'lat'].values
 
        ed['xc'] =  ed[['x1', 'x2']].mean(axis=1) #mean of the edge index
        ed['yc'] =  ed[['y1', 'y2']].mean(axis=1)

        ## min bottom index
        ed['kbs1'] = grd.loc[ed['node1'] - 1,'kbp00'].values
        ed['kbs2'] = grd.loc[ed['node2'] - 1,'kbp00'].values

        ed['kbs'] = ed[['kbs1', 'kbs2']].min(axis=1)
        
        

        xsides = xr.Dataset({
                         u'SCHISM_hgrid_edge_nodes' : ([u'nSCHISM_hgrid_edge', u'two'], [[x-1,y-1] for [x,y] in sides]), # index from 0
                         u'SCHISM_hgrid_edge_x' : ([u'nSCHISM_hgrid_edge'], ed['xc'].values),
                         u'SCHISM_hgrid_edge_y' : ([u'nSCHISM_hgrid_edge'], ed['yc'].values ),
                         u'edge_bottom_index' : ([u'nSCHISM_hgrid_edge'], ed.kbs.values)})                 

        logger.info('done with side based variables \n')

        # General properties
        
        header2 = self.misc['header'].apply(pd.to_numeric)
        nlist = ['start_year','start_month','start_day','start_hour','utc_start','dtout','nspool','nvrt','kz','ics']
        ddf = pd.DataFrame(header2).T
        ddf[nlist] = ddf[nlist].astype(int)
        header2 = ddf
        sigmas = [x for x  in header2.columns if 'sigma' in x]
        sigms = header2.loc[:,sigmas].values.flatten() # get sigmas
        iwet_dry = 0  # defined by the user
        ihgrid_id = -2147483647 # defined by user - 0,dummy_dim,ihgrid_id
        one = xr.Dataset({'dry_value_flag': (('one'), [iwet_dry]),'SCHISM_hgrid': (('one'),[ihgrid_id]) })

        #compute cs
        klev = np.arange(header2.kz.values[0],header2.nvrt.values[0]+1)
        k = klev-header2.kz.values


        cs=np.zeros(k)

        cs=(1-header2.theta_b.values)*np.sinh(header2.theta_f.values*sigms[k])/np.sinh(header2.theta_f.values)+ \
            header2.theta_b.values*(np.tanh(header2.theta_f.values*(sigms[k]+0.5))-np.tanh(header2.theta_f.values*0.5))/2/np.tanh(header2.theta_f.values*0.5)

        Cs = xr.Dataset({'Cs': (('sigma'), cs)}, coords={'sigma':sigms})

        header_list = ['ics','h0','h_c','theta_b','theta_f','h_s']
        gen = header2[header_list].to_xarray()
        gen = gen.rename({'index':'one'})

        
        #merge
        gen = xr.merge([gen,Cs,one])

        gen = gen.rename({'ics':'coordinate_system_flag','h0':'minimum_depth','h_c':'sigma_h_c','theta_b':'sigma_theta_b','theta_f':'sigma_theta_f','h_s':'sigma_maxdepth'})
        
        gen = gen.drop_vars('one')
        
        #set timestamp
        date = header2.loc[:,['start_year','start_month','start_day','start_hour','utc_start']]
        date = date.astype(int)
        date.columns=['year','month','day','hour','utc'] # rename the columns
        #set the start timestamp
        sdate = pd.Timestamp(year=date.year.values[0], month=date.month.values[0], day=date.day.values[0], hour=date.hour.values[0], tz=date.utc.values[0])

        logger.info('done with generic variables \n')
        
        # Read Netcdf output files
        hfiles = glob.glob(path+'outputs/schout_*_*.nc')
        
        irange = [int(x.split('_')[-1].split('.')[0]) for x in hfiles]
        irange = np.unique(irange)
        
        
        self.read_vgrid() # read grid attributes
        
        
#        total_xdat = []
        for val in irange:
            hfiles=glob.glob(path+'outputs/schout_*_{}.nc'.format(val))
            hfiles.sort()
            
            times = xr.open_dataset(hfiles[0]).time
            times = pd.to_datetime(times.values, unit='s',
                       origin=sdate.tz_convert(None))
            
            if times.size == 0 : continue
            
            idat = self.tcombine(hfiles,sdate, times)
                    
            #MERGE
                
            xc = xr.merge([idat,gen,xnodes,xelems,xsides])

            #Choose attrs
            if header2.ics.values == 1:
                lat_coord_standard_name = 'projection_y_coordinate'
                lon_coord_standard_name = 'projection_x_coordinate'
                x_units = 'm'
                y_units = 'm'
                lat_str_len = 23
                lon_str_len = 23
            else:
                lat_coord_standard_name = 'latitude'
                lon_coord_standard_name = 'longitude'
                x_units = 'degrees_east'
                y_units = 'degrees_north'
                lat_str_len = 8
                lon_str_len = 9
            
            #set Attrs
            xc.SCHISM_hgrid_node_x.attrs = {'long_name' : 'node x-coordinate', 'standard_name' : lon_coord_standard_name , 'units' : x_units, 'mesh' : 'SCHISM_hgrid'}

            xc.SCHISM_hgrid_node_y.attrs = {'long_name' : 'node y-coordinate', 'standard_name' : lat_coord_standard_name , 'units' : y_units, 'mesh' : 'SCHISM_hgrid'}

            xc.depth.attrs = {'long_name' : 'Bathymetry', 'units' : 'meters', 'positive' : 'down', 'mesh' : 'SCHISM_hgrid', 'location' : 'node'}

            xc.sigma_h_c.attrs = {'long_name' : 'ocean_s_coordinate h_c constant', 'units' : 'meters', 'positive' : 'down'}

            xc.sigma_theta_b.attrs = {'long_name' : 'ocean_s_coordinate theta_b constant'}

            xc.sigma_theta_f.attrs = {'long_name' : 'ocean_s_coordinate theta_f constant'}

            xc.sigma_maxdepth.attrs = {'long_name' : 'ocean_s_coordinate maximum depth cutoff (mixed s over z boundary)', 'units' : 'meters', 'positive' : 'down'}

            xc.Cs.attrs = {'long_name' : 'Function C(s) at whole levels', 'positive' : 'up' }

            xc.dry_value_flag.attrs = {'values' : '0: use last-wet value; 1: use junk'}

            xc.SCHISM_hgrid_face_nodes.attrs = {'long_name' : 'Horizontal Element Table', 'cf_role' : 'face_node_connectivity' , 'start_index' : 0}

            xc.SCHISM_hgrid_edge_nodes.attrs = {'long_name' : 'Map every edge to the two nodes that it connects', 'cf_role' : 'edge_node_connectivity' , 'start_index' : 0}

            xc.SCHISM_hgrid_edge_x.attrs = {'long_name' : 'x_coordinate of 2D mesh edge' , 'standard_name' : lon_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xc.SCHISM_hgrid_edge_y.attrs = {'long_name' : 'y_coordinate of 2D mesh edge' , 'standard_name' : lat_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xc.SCHISM_hgrid_face_x.attrs = {'long_name' : 'x_coordinate of 2D mesh face' , 'standard_name' : lon_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xc.SCHISM_hgrid_face_y.attrs = {'long_name' : 'y_coordinate of 2D mesh face' , 'standard_name' : lat_coord_standard_name, 'units' : 'm', 'mesh' : 'SCHISM_hgrid'}

            xc.SCHISM_hgrid.attrs = {'long_name' : 'Topology data of 2d unstructured mesh',
                                       'topology_dimension' : 2,
                                       'cf_role' : 'mesh_topology',
                                       'node_coordinates' : 'SCHISM_hgrid_node_x SCHISM_hgrid_node_y',
                                       'face_node_connectivity' : 'SCHISM_hgrid_face_nodes',
                                       'edge_coordinates' : 'SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y',
                                       'face_coordinates' : 'SCHISM_hgrid_face_x SCHISM_hgrid_face_y',
                                       'edge_node_connectivity' : 'SCHISM_hgrid_edge_nodes'
                                      }

            xc.node_bottom_index.attrs = {'long_name' : 'bottom level index at each node' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'node',
                'start_index' : 0}

            xc.ele_bottom_index.attrs = {'long_name' : 'bottom level index at each element' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'elem',
                'start_index' : 0}

            xc.edge_bottom_index.attrs = {'long_name' : 'bottom level index at each edge' , 'units' : 'non-dimensional', 'mesh' : 'SCHISM_hgrid', 'location' : 'edge',
                'start_index' : 0}
            
        
            base_date = ' '.join([str(x) for x in date.T.values.flatten()])
            xc.time.attrs = {'long_name': 'Time', 'base_date' : base_date , 'standard_name' : 'time' }
                
            xc.sigma.attrs ={'long_name' : 'S coordinates at whole levels',
                         'units' : '1',
                         'standard_name' : 'ocean_s_coordinate',
                         'positive' : 'up',
                         'h_s' : self.misc['hs'],
                         'h_c' : self.misc['h_c'],
                         'theta_b' : self.misc['theta_b'],
                         'theta_f' : self.misc['theta_f'],
                         'formula_terms':
                         's: sigma eta: elev depth: depth a: sigma_theta_f b: sigma_theta_b depth_c: sigma_h_c'}
            
            # Dataset Attrs

            xc.attrs = {'Conventions': 'CF-1.0, UGRID-1.0', 'title': 'SCHISM Model output', 'source': 'SCHISM model output version v10', 'references': 'http://ccrm.vims.edu/schismweb/',
                         'history': 'created by pyPoseidon', 'comment': 'SCHISM Model output', 'type': 'SCHISM Model output', 'VisIT_plugin': 'https://schism.water.ca.gov/library/-/document_library/view/3476283' }
        
        
            xc.to_netcdf(path+'outputs/schout_{}.nc'.format(val))
        
        
        
        logger.info('done with output netCDF files \n')
                


    def set_obs(self,**kwargs):
        
        path = get_value(self,kwargs,'rpath','./') 
        nspool_sta = get_value(self,kwargs,'nspool_sta',1)
        tg_database = get_value(self,kwargs,'tide_gauges',None) # TODO
        coastal_monitoring = get_value(self,kwargs,'coastal_monitoring',False)
        flags = get_value(self,kwargs,'station_flags',[1]+[0]*8) 
        
        station_flag = pd.DataFrame({'elev':flags[0], 'air_pressure':flags[1], 'windx':flags[2], 'windy':flags[3], 'T':flags[4], 'S':flags[5], 'u':flags[6], 'v':flags[7], 'w':flags[8]}, index =[0])
        
        
        z = self.__dict__.copy()
        ## FOR TIDE GAUGE MONITORING
        tg = obs.obs(**z)
        logger.info('get in-situ measurements locations \n')
        
        gpoints = np.array(list(zip(self.grid.Dataset.SCHISM_hgrid_node_x.values,self.grid.Dataset.SCHISM_hgrid_node_y.values)))
        
        stations = []
        grid_index = []
        for l in range(tg.locations.shape[0]):
            plat,plon = tg.locations.loc[l,['lat','lon']]
            cp = closest_node([plon,plat],gpoints)
            grid_index.append(list(zip(self.grid.Dataset.SCHISM_hgrid_node_x.values,self.grid.Dataset.SCHISM_hgrid_node_y.values)).index(tuple(cp)))
            stations.append(cp)
            
        #to df
        stations = pd.DataFrame(stations,columns=['SCHISM_hgrid_node_x','SCHISM_hgrid_node_y'])
        stations['z']=0
        stations.index += 1
        stations['gindex'] = grid_index
        
        if coastal_monitoring :
        ## FOR COASTAL MONITORING
            logger.info('set land boundaries as observation points \n')
        
            #get land boundaries
            coasts = [key for key in self.grid.Dataset.variables if 'land_boundary' in key]
            # get index
            bnodes=[]
            for i in range(len(coasts)):
                bnodes = bnodes + list(self.grid.Dataset[coasts[i]].dropna('index').astype(int).values)
            bnodes = np.unique(bnodes) # keep unique values
            # get x,y
            xx = self.grid.Dataset.SCHISM_hgrid_node_x.to_dataframe().loc[bnodes]
            yy = self.grid.Dataset.SCHISM_hgrid_node_y.to_dataframe().loc[bnodes]
            #put them in dataframe
            coastal_stations = pd.concat([xx, yy], axis=1, sort=False)
            coastal_stations['z'] = 0
            coastal_stations.reset_index(inplace=True, drop=True)
            coastal_stations.index += 1
            coastal_stations['gindex'] = bnodes
            coastal_stations.head()
            #append
            stations = pd.concat([stations,coastal_stations])
            stations.reset_index(inplace=True,drop=True)
            stations.index += 1            

        
        stations['gindex'] = stations['gindex'].astype(int)
        
        #modify config paramater
        self.params.loc['iout_sta'] = 1
        self.params.loc['nspool_sta'] = nspool_sta
        self.params.to_csv(path + 'param.in', header=None, sep='=') 
                
        self.stations = stations
        
        logger.info('write out stations.in file \n')
        
        #output to file
        with open(path + 'station.in', 'w') as f:
            station_flag.to_csv(f,header=None,index=False,sep=' ')
            f.write('{}\n'.format(stations.shape[0]))
            stations.loc[:,['SCHISM_hgrid_node_x',	'SCHISM_hgrid_node_y',	'z']].to_csv(f, header=None, sep=' ')
    
    def get_obs(self,**kwargs):
        
        # locate the station files
        sfiles = glob.glob(self.rpath + 'outputs/staout_*')
        sfiles.sort()

        try:
            # get the station flags
            flags = pd.read_csv(self.rpath + 'station.in', header=None,nrows=1,delim_whitespace=True).T
            flags.columns = ['flag']
            flags['variable'] = ['elev', 'air_pressure', 'windx', 'windy', 'T', 'S', 'u', 'v', 'w']
        
            vals = flags[flags.values==1]# get the active ones
        except OSError as e:
            if e.errno == errno.EEXIST: logger.error('No station.in file present')
            return
            

        dfs=[]
        for idx in vals.index:
            obs = np.loadtxt(sfiles[idx])
            df = pd.DataFrame(obs)
            df = df.set_index(0)
            df.index.name = 'time'
            df.columns.name = vals.loc[idx,'variable']
            df.index = pd.to_datetime(self.start_date) + pd.to_timedelta(df.index,unit='S')
            pindex = pd.MultiIndex.from_product([df.T.columns,df.T.index])

            dfs.append(pd.DataFrame(df.values.flatten(),index = pindex, columns=[vals.loc[idx,'variable']]).to_xarray().rename({'level_0':'time','level_1':'point'}))        
        
        
        self.time_series = xr.combine_by_coords(dfs)