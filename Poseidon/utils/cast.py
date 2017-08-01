import numpy as np
import datetime
import sys
import os
from shutil import copy2
import logging
import glob
import pickle
from Poseidon.model import *

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
        
                     
        files=['config_d_hydro.xml','*.grd','*.enc','*.obs','*.dep', '*.bnd', '*.bca','run_flow2d3d.sh','info.pkl']
        
                
        prev=self.folders[0]
        fpath = self.path+'/{}/'.format(prev)
        
        
        cf = [glob.glob(self.path+prev+'/'+e) for e in files]
        cfiles = [item.split('/')[-1] for sublist in cf for item in sublist]
                    
 
        for date,folder,meteo in zip(self.dates[1:],self.folders[1:],self.meteo_files[1:]):
            
            ppath = self.path+'/{}/'.format(prev)
            if not os.path.exists(ppath):
                sys.stdout.write('Initial folder not present {}\n'.format(ppath)) 
                sys.exit(1)
            
            prev = folder    
            # create the folder/run path

            rpath=self.path+'/{}/'.format(folder)   

            if not os.path.exists(rpath):
                os.makedirs(rpath)

            # load model
            with open(ppath+'info.pkl', 'r') as f:
                          info=pickle.load(f)
            
#            for attr, value in self.info.iteritems():
#                setattr(m, attr, value)
            m=model(**info)
                               
            #update the properties 
            m.impl.date = date
            m.impl.model['date'] = date
            m.impl.mpath=meteo 
            m.impl.model['mpath'] = meteo
            m.impl.rpath=rpath 
            m.impl.model['rpath'] = rpath
            
            
            # copy/link necessary files

            for filename in cfiles:
        #        copy2(ppath+filename,rpath+filename)
                os.symlink(fpath+filename,rpath+filename)
            
            copy2(ppath+m.impl.tag+'.mdf',rpath) #copy the mdf file
                
            # copy restart file

            inresfile='tri-rst.'+m.impl.tag+'.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

            outresfile='restart.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

          #  copy2(ppath+inresfile,rpath+'tri-rst.'+outresfile)
            os.symlink(ppath+inresfile,rpath+'tri-rst.'+outresfile)

            #get new meteo 

            sys.stdout.write('process meteo\n')
            sys.stdout.flush()

            
            check=[os.path.exists(rpath+f) for f in ['u.amu','v.amv','p.amp']]   
            if np.any(check)==False :
               
                m.force()
                m.uvp()  #write u,v,p files 
        
            else:
                sys.stdout.write('meteo files present\n')
            
            
            # modify mdf file    
            inp, order = mdf.read(rpath+m.impl.tag+'.mdf')    
            
            
            # adjust iteration date
            tstart = (date.hour+m.impl.ft1)*60
            tend = tstart + (m.impl.ft2)*60
            
            inp['Itdate']=datetime.datetime.strftime(date,'%Y-%m-%d')
            inp['Tstart']=[tstart]
            inp['Tstop']=[tend]
            inp['Flmap']  = [tstart,60,tend]
            inp['Flhis']  = [tstart,inp['Dt'][0],tend]


            if 'Restid' not in order : order.append('Restid')
              # adjust restart file   
            inp['Restid']=outresfile

            # update mdf
            mdf.write(inp, rpath+m.impl.tag+'.mdf',selection=order)
                                  
            # run case
            sys.stdout.write('executing\n')
            sys.stdout.flush()
         
            os.chdir(rpath)
            #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
            m.run()

            m.save()
            
            logging.info('done for date :'+datetime.datetime.strftime(date,'%Y%m%d.%H'))
            
