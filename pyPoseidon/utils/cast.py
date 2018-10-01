"""
Simulation management module

"""
# Copyright 2018 European Union
# This file is part of [software name], a software written by [author's name] ([JRC unit])
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
# See the Licence for the specific language governing permissions and limitations under the Licence. 

import numpy as np
import datetime
import sys
import os, errno
from shutil import copy2
import logging
import glob
import pickle
import pyPoseidon.model as pm
from pyPoseidon.utils.get_value import get_value
import pandas as pd
from pyPoseidon.utils import data

class cast:
    impl=None
    def __init__(self,**kwargs):
        model = kwargs.get('solver', None)
        if model == 'd3d' :
            self.impl = dcast(**kwargs)
        elif model == 'schism' :
            self.impl = scast(**kwargs)
        
    def run(self,**kwargs):
        self.impl.run(**kwargs)
    

class dcast(cast):
    
    def __init__(self,**kwargs):
               
        for attr, value in kwargs.iteritems():
                setattr(self, attr, value)
                
        logging.basicConfig(filename=self.path+self.case+'.log',level=logging.INFO)            
                   
    def run(self,**kwargs):
        
                      
        files=[self.tag+'_hydro.xml',self.tag+'.enc',self.tag+'.obs', self.tag+'.bnd', self.tag+'.bca','run_flow2d3d.sh']
        files_sym=[self.tag+'.grd',self.tag+'.dep']
        
                
        prev=self.folders[0]
        fpath = self.path+'/{}/'.format(prev)
        
        cf = [glob.glob(self.path+prev+'/'+e) for e in files]
        cfiles = [item.split('/')[-1] for sublist in cf for item in sublist]
                    
        for date,folder,meteo,time_frame in zip(self.dates[1:],self.folders[1:],self.meteo_files[1:],self.time_frames[1:]):
            
            ppath = self.path+'/{}/'.format(prev)
            if not os.path.exists(ppath):
                sys.stdout.write('Initial folder not present {}\n'.format(ppath)) 
                sys.exit(1)
            
            prev = folder    
            # create the folder/run path

            rpath=self.path+'/{}/'.format(folder)   

            if not os.path.exists(rpath):
                os.makedirs(rpath)

            copy2(ppath+self.tag+'_info.pkl',rpath) #copy the info file

            # load model
            with open(rpath+self.tag+'_info.pkl', 'r') as f:
                          info=pickle.load(f)
            
            args = set(kwargs.keys()).intersection(info.keys()) # modify dic with kwargs
            for attr in list(args):
                info[attr] = kwargs[attr]
            
            #update the properties   
            info['date'] = date
            info['start_date'] = date
            info['time_frame'] = time_frame
            info['mpaths'] = meteo
            info['rpath'] = rpath
            if self.rstep:
                info['rstep'] = self.rstep
            
#            for attr, value in self.iteritems():
#                setattr(info, attr, value)
            m=pm.model(**info)
                                                         
            # copy/link necessary files

            for filename in cfiles:
                 copy2(ppath+filename,rpath+filename)
        #     if os.path.exists(rpath+filename)==False: 
        #        os.symlink(fpath+filename,rpath+filename)
        
        
            #symlink the big files
            for filename in files_sym:
                ipath = glob.glob(self.path+self.folders[0]+'/'+filename)[0]
                try:
                    os.symlink(ipath,rpath+filename)
                except OSError, e:
                  if e.errno == errno.EEXIST:
                      sys.stdout.write('Restart link present\n')
                      sys.stdout.write('overwriting\n')
                      os.remove(rpath+filename)
                      os.symlink(ipath,rpath+filename)
            
            copy2(ppath+m.impl.tag+'.mdf',rpath) #copy the mdf file
                
            # copy restart file

            inresfile='tri-rst.'+m.impl.tag+'.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

            outresfile='restart.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

          #  copy2(ppath+inresfile,rpath+'tri-rst.'+outresfile)
            try:
              os.symlink(ppath+inresfile,rpath+'tri-rst.'+outresfile)
            except OSError, e:
              if e.errno == errno.EEXIST:
                  sys.stdout.write('Restart link present\n')
                  sys.stdout.write('overwriting\n')
                  os.remove(rpath+'tri-rst.'+outresfile)
                  os.symlink(ppath+inresfile,rpath+'tri-rst.'+outresfile)
              else:
                  raise e            

            #get new meteo 

            sys.stdout.write('process meteo\n')
            sys.stdout.flush()

            flag = get_value(self,kwargs,'update',[])
            
            check=[os.path.exists(rpath+f) for f in ['u.amu','v.amv','p.amp']]

            if (np.any(check)==False) or ('meteo' in flag):
               
                m.force()
                m.impl.to_force(m.impl.meteo.impl.uvp,vars=['msl','u10','v10'],rpath=rpath)  #write u,v,p files 
        
            else:
                sys.stdout.write('meteo files present\n')
            
            # modify mdf file
            m.config(config_file = ppath+m.impl.tag+'.mdf', config={'Restid':outresfile}, output=True)
                                              
            # run case
            sys.stdout.write('executing\n')
            sys.stdout.flush()
         
            os.chdir(rpath)
            #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
            m.run()

            m.save()
            
            #cleanup
            os.remove(rpath+'tri-rst.'+outresfile)
            
            # save compiled nc file
            
            out = data(**{'solver':m.impl.solver,'rpath':rpath,'savenc':True})
            
            logging.info('done for date :'+datetime.datetime.strftime(date,'%Y%m%d.%H'))
            
