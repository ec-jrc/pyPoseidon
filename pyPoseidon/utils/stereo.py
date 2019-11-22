import numpy as np


def to_stereo(lon,lat, R=1):

    # to 3D
    kx=np.cos(lat/180*np.pi)*np.cos(lon/180*np.pi)*R
    ky=np.cos(lat/180*np.pi)*np.sin(lon/180*np.pi)*R
    kz=np.sin(lat/180*np.pi)*R 
    
    # to 2D
    u = 2*R*kx/(R+kz)
    v = 2*R*ky/(R+kz)
    
    return u,v
    