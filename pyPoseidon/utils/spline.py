import pandas as pd
import geopandas as gp
import numpy as np
import shapely

def spline(df):

    pol = gp.GeoDataFrame(geometry=[shapely.geometry.LineString(df.values)])

    nls = len(pol.geometry[0].coords.xy[0])
    nls

    ips = 2*nls
    ips

    gpf= gp.GeoDataFrame(geometry=[pol.geometry[0].interpolate(i/float(ips-1), normalized=True) for i in range(ips)])

    di=pd.DataFrame(np.array([np.array(p.geometry.xy[:]).flatten() for (l,p) in gpf.iterrows()]),columns=['lon','lat'])


    return di