"""
Main model module of pyPoseidon. It controls the creation, output & execution of a complete simulation based on different hydrological models

Currently supported : DELFT3D , SCHISM
               
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

#local modules
from .d3d import *
from .schism import *
from pyPoseidon.bnd import *
import pyPoseidon
import pyPoseidon.grid as pgrid
import pyPoseidon.meteo as pmeteo
import pyPoseidon.dem as pdem
from pyPoseidon.utils.get_value import get_value
import logging
import colorlog

LOGGING = {
    'version': 1,
    'disable_existing_loggers': True,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
        'color': {
            '()': 'colorlog.ColoredFormatter',
        }
    },
    'handlers': {
        'file': {
            'level':'DEBUG',
            'class':'logging.FileHandler',
            'formatter': 'verbose',
            'filename': 'model.log'
        }
    ,
        'console':{
            'level':'DEBUG',
            'class':'colorlog.StreamHandler',
            'formatter': 'color'
        }
    },
    'loggers': {
        'pyPoseidon': {
            'handlers':['file','console'],
            'propagate': True,
            'level':'DEBUG',
        }
    }
}

logging.config.dictConfig(LOGGING)

#retrieve the module path
#DATA_PATH = pkg_resources.resource_filename('pyPoseidon', 'misc')
DATA_PATH = os.path.dirname(pyPoseidon.__file__)+'/misc/'    

# strings to be used 
le=['A','B']

nm = ['Z', 'A']


def model(solver=None, atm=True, tide=False, **kwargs):
    """
    Construct a hydrodynamic model based on different solvers.

    :param solver: Type of model, e.g. 'd3d' or 'schism'
    :param kwargs: additional arguments to pass along
    :type solver: String
    :type kwargs: dict
    :return: A model
    :rtype: composite object


    Example:

    model = pyPoseidon.model(solver='d3d)

    """
    
    kwargs.update({'atm':atm,'tide':tide})
    if solver == 'd3d':
        return d3d(**kwargs)
    elif solver == 'schism':
        return schism(**kwargs)            
        
def read_model(filename,**kwargs):
            
    end = filename.split('.')[-1]

    if end in ['txt','json'] :
        with open(filename, 'rb') as f:
            info = pd.read_json(f,lines=True).T
            info[info.isnull().values] = None
            info = info.to_dict()[0]
    else:
        logger.error('Model file should be .txt or .json')
        sys.exit(0)
    
    return model(**info)



        

