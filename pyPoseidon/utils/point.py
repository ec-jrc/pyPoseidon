"""
Point analysis module

"""
# Copyright 2018 European Union
# This file is part of pyPoseidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.


import numpy as np
from pyPoseidon.utils.pplot import *
from pyPoseidon.utils.obs import obs
from pyPoseidon.grid import *
import pyPoseidon.model as pmodel
import pyresample
import pandas as pd
import datetime
from pyPoseidon.utils.get_value import get_value
import holoviews as hv
import geoviews as gv
from cartopy import crs
import xarray as xr
import glob
import sys
import subprocess
import scipy.interpolate


class point:

    def __init__(self, solver = None, lon=None, lat = None, dataset = None, **kwargs):

        self.lon = lon
        self.lat = lat
        self.data = dataset

        if solver == 'd3d':
            self.time_series = self.tseries(**kwargs)
        elif solver == 'schism':
            self.time_series = self.stseries(**kwargs)
        else:
            logger.error('solver is not defined, exiting \n')
            sys.exit(1)


    def tseries(self,**kwargs):

        var = kwargs.get('var', None)
        method = kwargs.get('method', 'nearest')

        plat=float(self.lat)
        plon=float(self.lon)

        X, Y = self.data.Dataset.XZ.values[1:-1,1:-1].T,self.data.Dataset.YZ.values[1:-1,1:-1].T
        xh = np.ma.masked_array(X, self.data.w) #mask land
        yh = np.ma.masked_array(Y, self.data.w)


        i=np.abs(X[0,:]-plon).argmin()
        j=np.abs(Y[:,0]-plat).argmin()

        xb, yb = xh[j-5:j+5,i-5:i+5],yh[j-5:j+5,i-5:i+5]

        vals = self.data.Dataset[var][:,j-5:j+5,i-5:i+5].values

        orig = pyresample.geometry.SwathDefinition(lons=xb,lats=yb) # create original swath grid

        targ = pyresample.geometry.SwathDefinition(lons=np.array([plon,plon]),lats=np.array([plat,plat])) #  point

        svals = []

        if method == 'nearest':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_nearest(orig,vals[k,:,:],targ,radius_of_influence=100000,fill_value=np.nan)
                svals.append(s[0])
        elif method == 'gauss':
            for k in range(vals.shape[0]):
                s = pyresample.kd_tree.resample_gauss(orig,vals[k,:,:],targ,radius_of_influence=100000,fill_value=np.nan,sigmas=25000)
                svals.append(s[0])


        pdata = pd.DataFrame({'time':self.data.Dataset[var].time, self.data.Dataset[var].name : svals})
        return pdata.set_index(['time'])


    def stseries(self,**kwargs):

        var = kwargs.get('var', None)
        method = kwargs.get('method', 'nearest')

        plat=float(self.lat)
        plon=float(self.lon)

        points = pd.concat([self.data.SCHISM_hgrid_node_x[:].to_dataframe(), self.data.SCHISM_hgrid_node_y[:].to_dataframe()], axis=1)
        values = self.data[var].values

        svals = []
        for k in range(self.data[var].time.size):
            svals.append( scipy.interpolate.griddata(points.values, values[k,:], (plon, plat), method=method) )

        pdata = pd.DataFrame({'time':self.data[var].time, self.data[var].name : svals})
        return pdata.set_index(['time'])

