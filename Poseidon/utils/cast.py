import numpy as np
import datetime
import sys
import os
from shutil import copy2
import logging
import glob
import pickle


class cast:
            
    def __init__(self,**kwargs):
        
        self.path = kwargs.get('path', None)
        self.case = kwargs.get('case', None)
    
        logging.basicConfig(filename=self.path+self.case+'.log',level=logging.INFO)
    
        self.sd = kwargs.get('start_date', None)
        self.fd = kwargs.get('end_date', None)
        self.dt = kwargs.get('time_interval', None)
        
        self.dates = kwargs.get('date_list', None)
        self.folders = kwargs.get('folders', None)
        self.meteo = kwargs.get('meteo_files', None)

    
    def d3d(self,**kwargs):
        
        
                
        files=['config_d_hydro.xml','*.mdf','*.grd','*.enc','*.obs','*.dep', '*.bnd', '*.bca','run_flow2d3d.sh']
        
        cf = [glob.glob(e) for e in files]
        cfiles = [item for sublist in cf for item in sublist]
        
        prev=self.folders[0]
        pdate=self.dates[0]
        
        with open(self.path+prev+'/info.pkl', 'r') as f:
            
              self.info=pickle.load(f)
        
        
        for date,folder,meteo in zip(self.dates[1:],self.folders[1:],self.meteo[1:]):
            
            ppath = self.path+'/{}/'.format(prev)
            if not os.path.exists(fpath):
                sys.stdout.write('Initial folder not present {}\n'.format(fpath)) 
                sys.exit(1)
                
            # create the folder/run path

            rpath=self.path+'/{}/'.format(folder)   

            if not os.path.exists(rpath):
                os.makedirs(rpath)

            # copy necessary files

            for filename in cfiles:
                copy2(ppath+filename,rpath+filename)

            # copy restart file

            inresfile='tri-rst.'+self.info.tag+'.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

            outresfile='restart.'+datetime.datetime.strftime(date,'%Y%m%d.%H%M%M')

            copy2(ppath+inresfile,rpath+'tri-rst.'+outresfile)

            #get new meteo 

            sys.stdout.write('process meteo\n')
            sys.stdout.flush()

            check=[os.path.exists(rpath+f) for f in ['u.amu','v.amv','p.amp']]   
            if np.any(check)==False :

                m = ecmwf() # initialize
                m.parse(path=meteo,**self.info)
                m.force(path=rpath)  #write u,v,p files 
        
            else:
                sys.stdout.write('meteo files present\n')
            
            
            # modify mdf file    
            inp, order = mdf.read(rpath+bname+'.mdf')    
            
            # adjust iteration date
            tstart=rundate.hour*60
            inp['Itdate']=datetime.datetime.strftime(date,'%Y-%m-%d')
            inp['Tstart']=[tstart]
            inp['Tstop']=[nt*60+tstart]
            inp['Flmap']  = [tstart,60,nt*60+tstart]
            inp['Flhis']  = [tstart,inp['Dt'][0],nt*60+tstart]


            if 'Restid' not in order : order.append('Restid')
              # adjust restart file   
            inp['Restid']=outresfile

            # update mdf
            mdf.write(inp, rpath+self.tag+'.mdf',selection=order)
            # run case
         
            os.chdir(rpath)
            #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
            os.system('./run_flow2d3d.sh')

            
            logging.info('done for date :'+datetime.datetime.strftime(idate,'%Y%m%d.%H'))
         