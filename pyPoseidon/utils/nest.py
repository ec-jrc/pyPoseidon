import numpy as np
import datetime
import sys
import os, errno
from shutil import copy2
import logging
import glob
import pickle
from pyPoseidon.model import *
from pyPoseidon.util import *


class nest:
    impl=None
    def __init__(self,**kwargs):
        model = kwargs.get('solver', None)
        if model == 'd3d' :
            self.impl = dnest(**kwargs)
        elif model == 'schism' :
            self.impl = snest(**kwargs)
        
    def run(self,**kwargs):
        self.impl.run(**kwargs)
    

class dnest(nest):
    
    def __init__(self,**kwargs):
               
        for attr, value in kwargs.iteritems():
                setattr(self, attr, value)
                
        logging.basicConfig(filename=self.path+self.case+'.log',level=logging.INFO)    
        
        
        parent = kwargs.get('parent', '.')        
        
        # load model
        with open(parent+'info.pkl', 'r') as f:
                      info=pickle.load(f)
   
        
        info['minlon'] = kwargs.get('minlon', None)        
        info['maxlon'] = kwargs.get('maxlon', None)        
        info['minlat'] = kwargs.get('minlat', None)        
        info['maxlat'] = kwargs.get('maxlat', None)    
        
            
        info['rpath'] = kwargs.get('rpath', info['rpath']+'./nested/')
        rpath = info['rpath'] #save for later
        
        info['resolution'] = kwargs.get('resolution', None)
        
        info['atm'] = False
        
        
        #create new case 
        nest = model(**info)    
        
        nest.set() #setup nested run 
        
        nest.output() #output to run folder      
            
        check=[os.path.exists(parent+f) for f in ['u.amu','v.amv','p.amp']]   
        if np.any(check)==False :
               
                nest.force()
                nest.uvp()  #write u,v,p files 
        
        else: #link u,v,p
            for filename in ['u.amu','v.amv','p.amp']:
                os.symlink(parent+filename,rpath+filename)
                
                
        # modify mdf file    
        inp, order = mdf.read(rpath+nest.impl.tag+'.mdf')    
                    
        # adjust variables
        
        #create the ini file
        
        if 'Filic' not in order: order.append('Filic')
        inp['Filic']=nest.impl.tag + '.ini'
        
        pdata = data([parent])
        
        s1 = pdata.get_data('S1',step=1)
        u1 = pdata.get_data('U1',step=1)
        v1 = pdata.get_data('V1',step=1)
        
        xz=pdata.get_data('XZ')
        yz=pdata.get_data('YZ')
        
        orig = pyresample.geometry.SwathDefinition(lons=xz,lats=yz) # original points
        targ = pyresample.geometry.SwathDefinition(lons=nest.impl.grid.x,lats=nest.impl.grid.y) # target grid
        
        s_ini = pyresample.kd_tree.resample_nearest(orig,s1,targ,radius_of_influence=100000,fill_value=0)
        u_ini = pyresample.kd_tree.resample_nearest(orig,u1,targ,radius_of_influence=100000,fill_value=0)
        v_ini = pyresample.kd_tree.resample_nearest(orig,v1,targ,radius_of_influence=100000,fill_value=0)
        
        
        # Write .ini file
        with open(rpath+nest.impl.tag+'.ini', 'w') as f:
           np.savetxt(f,s_ini)
           np.savetxt(f,u_ini)
           np.savetxt(f,v_ini)
        
        
         #create the bc file
        
        if 'Filbnd' not in order: order.append('Filbnd')
        if 'Filana' not in order: order.append('Filana')
        
        inp['Filbnd']=nest.impl.tag+'.bnd'
        inp['Filana']=nest.impl.tag+'.bca'

        
#        bca = 


        inp['Restid']='##' #no restart file

        # update mdf
        mdf.write(inp, rpath+nest.impl.tag+'.mdf',selection=order)
                                  
        # run case
        sys.stdout.write('executing\n')
        sys.stdout.flush()
         
        os.chdir(rpath)
        #subprocess.call(rpath+'run_flow2d3d.sh',shell=True)
        nest.run()

        nest.save()
            
        logging.info('nesting run done for date :'+datetime.datetime.strftime(date,'%Y%m%d.%H'))
            

#    def get_ini():
    
    
#    def get_bca():
        
        