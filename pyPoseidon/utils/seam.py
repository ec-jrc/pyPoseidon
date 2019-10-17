import pandas as pd
import geopandas as gp
import numpy as np
import shapely
from shapely.ops import triangulate


def to_2d(x,y,z,tri3,**kwargs):
            
    if z==None : z = 1.
    
    gr = pd.DataFrame({'lon':x, 'lat':y, 'z':z})
    tf = gr[(gr.lon<-90) | (gr.lon>90.) | (gr.lat > 80)]
    
    elems = pd.DataFrame(tri3, columns=['a','b','c'])
    bes = elems[(elems.a.isin(tf.index.values)) & (elems.b.isin(tf.index.values)) & (elems.c.isin(tf.index.values))]
    
    qq = bes.copy()
    qq['ap'] = qq.apply(lambda x : tf.loc[x.a,['lon','lat']].values,axis=1)
    qq['bp'] = qq.apply(lambda x : tf.loc[x.b,['lon','lat']].values,axis=1)
    qq['cp'] = qq.apply(lambda x : tf.loc[x.c,['lon','lat']].values,axis=1)
    qq['geometry'] =  qq.apply(
             lambda x : shapely.geometry.Polygon([x.ap,x.bp,x.cp]),axis=1)
    
    ge=gp.GeoDataFrame(qq.geometry)
    gemask = ge.geometry.bounds.diff(axis=1,periods=2).maxx>300 # Avoid seam for the 2D graph
    
    mels = qq[gemask]
    
    
    adv = 0 
    p1=[]
    p2 = []
    for indx,vals in mels.iterrows():
    #    print(indx)
        pol = vals.geometry
        lons = [x for (x,y) in list(pol.exterior.coords)]
        incl = np.array(lons[:-1]).mean()
        if incl < 0.:
            lons = [x - 360. if x > 0. else x for x in pol.boundary.xy[0]]
        else:
            lons = [x + 360. if x < 0. else x for x in pol.boundary.xy[0]]
    
        nns = list(zip(lons, list(pol.boundary.xy[1]))) #points of recasted element

        npol = shapely.geometry.LineString(nns) # make a geometrical object (triangle) out of the nodes above 

        # create a meridian line
        if incl < 0.:
            l = shapely.geometry.LineString([[-180.,-90.],[-180.,90.]])
        else:
            l = shapely.geometry.LineString([[180.,-90.],[180.,90.]])
    
        cp = npol.intersection(l) # find the intersection with the element

        try:
            cpp = [(x.coords[0][0],x.coords[0][1]) for x in cp] # get cross nodes
        except:
            cpp = []

        de = pd.DataFrame(cpp + nns[:-1], columns = ['lon','lat'])

        de = de.drop_duplicates().reset_index(drop=True)
   
        points = shapely.geometry.MultiPoint(nns[:-1]+cpp)

        triangles = triangulate(points)

        tes = []
        for triangle in triangles:
            blon = [ x for (x,y) in list(triangle.exterior.coords)]
            blat = [ y for (x,y) in list(triangle.exterior.coords)]
            tes.append(de[(de['lon'].isin(blon)) & (de['lat'].isin(blat))].index.values)
    
        nels = pd.DataFrame(tes, columns=['a','b','c'])

        if cpp : 
            de = de.append(de.loc[:1],ignore_index=True) # replicate the meridian cross points
            de.lon[-2:] *= -1
    

        if de[de.lon > 180.].size>0:
            ids = de[de.lon > 180.].index.values # problematic nodes
            de.loc[de.lon > 180., 'lon'] -= 360.
        elif de[de.lon < -180.].size>0:
            ids = de[de.lon < -180.].index.values # problematic nodes
            de.loc[de.lon < -180., 'lon'] += 360.
    
        p1.append(de)


        if cpp :
            des = nels.loc[(nels.a.isin(ids)) | (nels.b.isin(ids)) | (nels.c.isin(ids)) ].copy()
            des.loc[des.a == 0, 'a'] = de.index[-2]
            des.loc[des.a == 1, 'a'] = de.index[-1]
            des.loc[des.b == 0, 'b'] = de.index[-2]
            des.loc[des.b == 1, 'b'] = de.index[-1]
            des.loc[des.c == 0, 'c'] = de.index[-2]
            des.loc[des.c == 1, 'c'] = de.index[-1]
            nels.loc[des.index] = des
    

        nels = nels + adv # reindex to global index

        p2.append(nels)
    
        adv = nels.values.max() + 1

    ng = pd.concat(p1)
    ng['z'] = 1 ## ADJUST
    
    ng.reset_index(inplace=True,drop=True)
    
    nge = pd.concat(p2)
    nge.reset_index(inplace=True,drop=True)
    
    si = gr.index[-1]
    ## drop the problematic elements
    ges = elems.drop(qq[gemask].index)
    ## append new nodes
    mes = gr.append(ng)
    ## Make new elements index global
    nges = nge + si + 1
    # append new elements
    ges = ges.append(nges)
    ges.reset_index(inplace=True, drop=True)
    mes.reset_index(inplace=True, drop=True)
    
    xx = mes.lon.values
    yy = mes.lat.values
    return xx,yy,ges.values
