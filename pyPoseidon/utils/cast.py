import numpy as np
import datetime
import sys
import os, errno
from shutil import copy2
import logging
import glob
import pickle
from pyPoseidon.model import *

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
        
                      
        files=[self.tag+'_hydro.xml',self.tag+'.grd',self.tag+'.enc',self.tag+'.obs',self.tag+'.dep', self.tag+'.bnd', self.tag+'.bca','run_flow2d3d.sh']
        
                
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

            copy2(ppath+self.tag+'_info.pkl',rpath) #copy the info file

            # load model
            with open(rpath+self.tag+'_info.pkl', 'r') as f:
                          info=pickle.load(f)
            
#            for attr, value in self.info.iteritems():
#                setattr(m, attr, value)
            m=model(**info)
            
            if type(meteo) is not list : meteo = [meteo] # make it a list so that it works for meteo 
            #update the properties 
            m.impl.date = date
            m.impl.model['date'] = date
            m.impl.mpaths=meteo 
            m.impl.model['mpaths'] = meteo
            m.impl.rpath=rpath 
            m.impl.model['rpath'] = rpath
            
            # copy/link necessary files

            for filename in cfiles:
                 copy2(ppath+filename,rpath+filename)
        #     if os.path.exists(rpath+filename)==False: 
        #        os.symlink(fpath+filename,rpath+filename)
            
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

            
            check=[os.path.exists(rpath+f) for f in ['u.amu','v.amv','p.amp']]   
            if np.any(check)==False :
               
                m.force()
                m.impl.to_force(m.impl.meteo.impl.uvp,vars=['msl','u10','v10'],rpath=rpath)  #write u,v,p files 
        
            else:
                sys.stdout.write('meteo files present\n')
            
            
            # modify mdf file    
            mdf = pd.read_csv(rpath+m.impl.tag+'.mdf',sep='=')    
            
            mdf = mdf.set_index(mdf.columns[0])
            
            mdfidx = mdf.index.str.strip() # store the stripped names
            
            # adjust iteration date
            tstart = (date.hour+m.impl.ft1)*60
            tend = tstart + (m.impl.ft2)*60
            
            dt = mdf.loc[mdf.index.str.contains('Dt')].values[0]
            
            mdf.loc[mdf.index.str.contains('Itdate')]='#{}#'.format(datetime.datetime.strftime(date,'%Y-%m-%d'))
            mdf.loc[mdf.index.str.contains('Tstart')]=tstart
            mdf.loc[mdf.index.str.contains('Tstop')]=tend
            mdf.loc[mdf.index.str.contains('Flmap')]='{} {} {}'.format(tstart,m.impl.step,tend)
            mdf.loc[mdf.index.str.contains('Flhis')]='{} {} {}'.format(tstart,dt,tend)

            if not 'Restid' in mdfidx: 
                mdf.reindex(mdf.index.values.tolist()+['Restid '])

            mdf.loc['Restid '] = outresfile # adjust restart file

            # update mdf
            mdf.to_csv(rpath+m.impl.tag+'.mdf',sep='=')
                                  
            # run case
            sys.stdout.write('executing\n')
            sys.stdout.flush()
         
            os.chdir(rpath)
            #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
            m.run()

            m.save()
            
            #cleanup
            os.remove(rpath+'tri-rst.'+outresfile)
            
            logging.info('done for date :'+datetime.datetime.strftime(date,'%Y%m%d.%H'))
            
