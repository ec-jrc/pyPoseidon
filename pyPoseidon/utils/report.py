# Create html for visualizing time series

import pyPoseidon
import pandas as pd
import geopandas as gp
import numpy as np
import branca
import shapely
import folium
from folium import plugins
import altair as alt
import os
import json
import logging

alt.data_transformers.disable_max_rows()

logger = logging.getLogger('pyPoseidon')

def source(tg,st):
    nearest = st.reindex(tg.index,method='nearest')
    

    pdata = pd.concat({'observation':tg,'simulation':nearest},names=['legend', 'date'])
    pdata = pd.DataFrame(pdata)
    pdata.columns=['height']
    pdata = pdata.reset_index()
    return pdata
    
    
def graph(source, title):
    
   
    # Create a selection that chooses the nearest point & selects based on x-value
    nearest = alt.selection(type='single', nearest=True, on='mouseover',
                            fields=['date'], empty='none')

    line = alt.Chart(source,width=550,height=300).mark_line(interpolate='basis').encode(
        x='date',
        y='height',
        color='legend'
    )

    # Transparent selectors across the chart. This is what tells us
    # the x-value of the cursor
    selectors = alt.Chart(source).mark_point().encode(
        x='date',
        opacity=alt.value(0),
    ).add_selection(
        nearest
    )

    # Draw points on the line, and highlight based on selection
    points = line.mark_point().encode(
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    )

    # Draw text labels near the points, and highlight based on selection
    text = line.mark_text(align='left', dx=5, dy=-5).encode(
        text=alt.condition(nearest, 'height', alt.value(' '))
    )

    # Draw a rule at the location of the selection
    rules = alt.Chart(source).mark_rule(color='gray').encode(
        x='date',
    ).transform_filter(
        nearest
    )

    # Put the five layers into a chart and bind the data
    graph = alt.layer(
        line, selectors, points, rules, text
    ).properties(
        width=550, height=300
    )
    
    graph.title = title
    return graph


def to_html(path='./',tag='schism'):

    logger.info('create html file for reporting')
    
    m = pyPoseidon.read_model(path + tag +'_model.json')
    m.get_data(online=True)

    d = m.data.Dataset

    vmin = d.elev.min().compute()
    vmax = d.elev.max().compute()

    x = d.SCHISM_hgrid_node_x.values
    y = d.SCHISM_hgrid_node_y.values
    tri = d.SCHISM_hgrid_face_nodes.values

    nodes = pd.DataFrame({'lon':x,'lat':y})
    tria = pd.DataFrame(tri ,columns=['a','b','c'])

    tria['ap'] = tria.apply(lambda x : nodes.loc[x.a,['lon','lat']].values,axis=1)
    tria['bp'] = tria.apply(lambda x : nodes.loc[x.b,['lon','lat']].values,axis=1)
    tria['cp'] = tria.apply(lambda x : nodes.loc[x.c,['lon','lat']].values,axis=1)

    tria['geometry'] =  tria.apply(
         lambda x : shapely.geometry.Polygon([x.ap,x.bp,x.cp]),axis=1)

   
    colormap = branca.colormap.LinearColormap(['green', 'yellow', 'red'], vmin=vmin.values, vmax=vmax.values)
    colormap.caption = 'Elevation'
      
    # geopandas
    gf_ = gp.GeoDataFrame(tria, crs={'init' :'epsg:4326'})    

    gf_ = gf_.drop(['a','b','c','ap','bp','cp'], axis=1)
    

    
    colormap = branca.colormap.LinearColormap(['green', 'yellow', 'red'], vmin=vmin.values, vmax=vmax.values)
    colormap.caption = 'Elevation'


    # MAXH
    mh = d.elev.max('time').values
        
    tria['value'] = tria.apply(lambda x : np.mean([mh[x.a],mh[x.b],mh[x.c]]),axis=1)

    gf_ = gp.GeoDataFrame(tria, crs={'init' :'epsg:4326'})    

    gf_ = gf_.drop(['a','b','c','ap','bp','cp'], axis=1)


    maxh = folium.GeoJson(
    gf_,
    style_function=lambda x: {
        "fillColor": colormap(x['properties']['value']),
        'weight':0,
        "fillOpacity": 1
    },
    )

    maxh.layer_name = 'Maximum Elevation'

    # GRID

    
    g = folium.GeoJson(
    gf_,
    style_function=lambda x: {
        "fillColor": "transparent",
        "color": "black",
        'weight':.2


        },
    )

    g.layer_name = 'Grid'

    #VALIDATION
    
    tgs = m.data.obs.locations.loc[m.data.obs.locations.Group=='TD UNESCO']
    
    points = tgs[['lat','lon']].values
    
    mlat,mlon = points.mean(axis=0)
    
    tg = m.data.obs.data
    tg.index.levels[1].name = 'date'
    tg.index.levels[0].name = 'location'
    tg = tg.drop(['Level m','Tide m'],axis=1)
    tg.columns=['height']

    toponames = tgs.Name.to_list()
    toponames = [x.strip(' ') for x in toponames]
    
    st = m.data.time_series.to_dataframe()
    st.index.levels[1].name = 'date'
    st.index.levels[0].name = 'location'
    st.columns=['height']
    st.index = st.index.set_levels(m.data.obs.locations.Name.str.strip(' '), level=0)
    st = st.loc[toponames]    

    mc = plugins.MarkerCluster()
    for name,latlon in zip(toponames,points):
        popup = folium.Popup()
        tg_ = tg.loc[name.strip(' ')]
        st_ = st.loc[name.strip(' ')]
        if tg_.shape[0] <= 1:
            tg_ = st_.copy()
            tg_.height = 0
        s = source(tg_,st_)
        gr = graph(s, name)
        vis = gr.to_json()
        folium.VegaLite(vis, height=350, width=750).add_to(popup)
        mc.add_child(folium.Marker(latlon, popup=popup))
    
    mc.layer_name = 'Validation'
    
    #FINAL MAP
    
    Map = folium.Map(location=(mlat, mlon), tiles='OpenStreetMap', date=m.date.__str__(), zoom_start=5)
    
    plugins.Fullscreen(
    position='topleft',
    title='Expand me',
    title_cancel='Exit me',
    force_separate_button=True
    ).add_to(Map)

    Map.add_child(maxh)
    Map.add_child(g)
    Map.add_child(mc)
    Map.add_child(folium.map.LayerControl())
    Map.add_child(colormap)


    Map.save(path + 'report.html')
    
    logger.info('... saved')
    
    return
    