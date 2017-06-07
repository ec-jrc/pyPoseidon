import numpy as np
import datetime
import sys
import os
from shutil import copy2
import logging



class cast:
    
    
    def __init__(**kwargs):
        
    self.path = kwargs.get('path', None)
    
    dt=(fd-sd).total_seconds()
    ndt=dt/(3600*12)
    ndt=np.int(ndt)+1

    for it in range(ndt): 
      idate=sd+datetime.timedelta(hours=12*it)
      logging.info(datetime.datetime.strftime(idate,'%Y%m%d.%H'))
      go(idate,path,att,TAT)
    
    
    logging.basicConfig(filename=self.path+case+'.log',level=logging.INFO)
