import numpy as np
#https://gitlab.onelab.info/gmsh/gmsh/-/blob/master/Geo/gmshSurface.cpp#L88

def to_3d(x,y,R=1):
    
    lon = np.array(x)
    lat = np.array(y)
    # to 3D
    kx=np.cos(lat/180*np.pi)*np.cos(lon/180*np.pi)*R
    ky=np.cos(lat/180*np.pi)*np.sin(lon/180*np.pi)*R
    kz=np.sin(lat/180*np.pi)*R
    
    return kx,ky,kz

def to_stereo(x,y, R=1):

    kx,ky,kz = to_3d(x,y,R)

    # to 2D in stereo
#    u = 2*R*kx/(R+kz)
#    v = 2*R*ky/(R+kz)
    u=-kx/(R+kz)
    v=-ky/(R+kz)

    return u,v

def stereo_to_3d(u, v, R=1):

    #to 3D
#    c=4*R**2/(u**2+v**2+4*R**2)
#    x=c*u
#    y=c*v
#    z=2*c*R-R

    rp2 = u**2+v**2
    x = -2 * R * u / (1 + rp2)
    y = -2 * R * v / (1 + rp2)
    z = R * (1 - rp2) / (1 + rp2)

    return x,y,z

def to_lat_lon(x, y, z=None):

    if z is None:
        x,y,z = to_3d(x,y)

    #to lat/lon
    rad = x**2+y**2+z**2
    rad = np.sqrt(rad)

    rlat = np.arcsin(z/rad)
    rlon = np.arctan2(y,x)

    rlat = rlat * 180 / np.pi
    rlon = rlon * 180 / np.pi

    return rlon,rlat



