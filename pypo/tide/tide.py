"""Setting the tidal boundary conditions"""
import numpy as np
from netCDF4 import Dataset
from scipy import interpolate


class tide:
    impl=None
    def __init__(self, **kwargs):
        tmodel = kwargs.get('tmodel', None)
        if tmodel == 'tpx07' :
            self.impl = tpx07(**kwargs)

class tpx07(tide):
    """Get tide components from TPX07
    """
    def __init__(self, **kwargs):
        self.properties = kwargs.get('properties', None)
        lons = kwargs.get('blons', None)
        lats  = kwargs.get('blats', None)

        filename=kwargs.get('tpath', None)

        dmed=Dataset(filename)

        lat=dmed['lat'][:]
        lon=dmed['lon'][:]

        tidal_c=dmed['tidal_constituents'][:]

        tidal_c=[''.join(k).upper().strip() for k in tidal_c]

        amp=dmed['tidal_amplitude_h']
        ph=dmed['tidal_phase_h']
        depth = dmed['depth']

        #adjust lon according to the grid window
        if lons.min()<= 0.1 :
          ii=np.abs(lon-180.).argmin()
          topo = lon[ii:]-360.
          lon = np.hstack([topo,lon[:ii]])
          topo = amp[ii:,:,:]
          amp = np.vstack([topo,amp[:ii,:,:]])
          topo = ph[ii:,:,:]
          ph = np.vstack([topo,ph[:ii,:,:]])
          topo = depth[ii:,:]
          dep = np.vstack([topo,depth[:ii,:]])


        phv=np.zeros([lons.shape[0],ph.shape[-1]])
        amv=np.zeros([lons.shape[0],amp.shape[-1]])

        k = 0
        for plon,plat in zip(lons,lats):

            i=np.abs(lon-np.float(plon)).argmin()
            j=np.abs(lat-np.float(plat)).argmin()

       #     if dep[i,j]==0.:


            xx = lon[i-1:i+2]
            yy = lat[j-1:j+2]

            for m in range(amp.shape[-1]):

              zz = amp[i-1:i+2,j-1:j+2,m]
              fa=interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)

              amv[k,m]= fa(plon,plat)

              zz = ph[i-1:i+2,j-1:j+2,m]
              fa=interpolate.RectBivariateSpline(xx,yy,zz,kx=2,ky=2)
              phv[k,m]= fa(plon,plat)

            k+=1


        self.phase=phv
        self.ampl=amv
        self.constituents = tidal_c


