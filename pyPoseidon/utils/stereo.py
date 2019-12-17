import numpy as np


def to_stereo(x,y, R=1):

    lon = np.array(x)
    lat = np.array(y)
    # to 3D
    kx=np.cos(lat/180*np.pi)*np.cos(lon/180*np.pi)*R
    ky=np.cos(lat/180*np.pi)*np.sin(lon/180*np.pi)*R
    kz=np.sin(lat/180*np.pi)*R 
    
    # to 2D
    u = 2*R*kx/(R+kz)
    v = 2*R*ky/(R+kz)
    
    return u,v
    
def to_lat_lon(u, v, R=1):
    
    #to 3D
    c=4*R**2/(u**2+v**2+4*R**2)
    
    x=c*u
    y=c*v
    z=2*c*R-R
    
    #to lat/lon
    rad = x**2+y**2+z**2
    rad = np.sqrt(rad)
    
    rlat = np.arcsin(z/rad)
    rlon = np.arctan2(y,x)
    
    rlat = rlat * 180 / np.pi
    rlon = rlon * 180 / np.pi
    
    return rlon,rlat
    
    
    