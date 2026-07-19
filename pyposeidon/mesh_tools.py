"""
Mesh adjustment functions

"""

# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import numpy as np
import pandas as pd
import xarray as xr
from tqdm.auto import tqdm

import pyposeidon.mesh as pmesh
from pyposeidon.dem_tools import fix, dem_range, resample


def split_dataframe(df, chunk_size=1000):
    chunks = list()
    num_chunks = len(df) // chunk_size + (1 if len(df) % chunk_size else 0)
    for i in range(num_chunks):
        chunks.append(df[i * chunk_size : (i + 1) * chunk_size].copy())
    return chunks


def split_boundaries(m, chunk=1000, **kwargs):

    d = m[["id", "bnode", "type"]].to_dataframe()  # get boundary info
    grouped = d.groupby("id")  # group

    # split
    dfn = []
    idx = 1001
    for name, group in grouped:
        if group.shape[0] > chunk:
            #        print(name)
            dfs = split_dataframe(group, chunk_size=chunk)
            for df in dfs:
                df["type"] = "land"
                df["id"] = idx
                idx += 1
            dfn.append(dfs)
            d.drop(d[d.id == name].index, inplace=True)

    dfn = [j for i in dfn for j in i]  # join

    new = pd.concat(dfn)  # concat

    d_new = pd.concat([d, new]).sort_values(["type", "id", "bnode"]).reset_index(drop=True)

    d_new.index.name = "bnodes"

    mnew = m.drop_vars(["bnode", "type", "id"]).drop_vars("bnodes")  # drop old

    mnew = xr.merge([mnew, d_new.to_xarray()])  # merge new

    return mnew


def populate(dem, perms, m, coastlines, buffer=0.0):

    for (i1, i2), (j1, j2) in tqdm(perms, total=len(perms)):

        lon1 = dem.longitude.data[i1:i2][0]
        lon2 = dem.longitude.data[i1:i2][-1]
        lat1 = dem.latitude.data[j1:j2][0]
        lat2 = dem.latitude.data[j1:j2][-1]

        # buffer lat/lon
        blon1 = lon1 - buffer
        blon2 = lon2 + buffer
        blat1 = lat1 - buffer
        blat2 = lat2 + buffer

        #    de = dem.sel(lon=slice(blon1,blon2)).sel(lat=slice(blat1,blat2))
        de = dem_range(dem, blon1, blon2, blat1, blat2)

        de, check, flag = fix(de, coastlines)

        # subset mesh
        indices_of_nodes_in_bbox = np.where(
            (m.SCHISM_hgrid_node_y >= lat1 - buffer / 2)
            & (m.SCHISM_hgrid_node_y <= lat2 + buffer / 2)
            & (m.SCHISM_hgrid_node_x >= lon1 - buffer / 2)
            & (m.SCHISM_hgrid_node_x <= lon2 + buffer / 2)
        )[0]

        bm = m.isel(nSCHISM_hgrid_node=indices_of_nodes_in_bbox)
        ids = np.argwhere(np.isnan(bm.depth.values)).flatten()
        #        grid_x, grid_y = bm.SCHISM_hgrid_node_x[ids], bm.SCHISM_hgrid_node_y[ids]
        grid_x, grid_y = bm.SCHISM_hgrid_node_x.data, bm.SCHISM_hgrid_node_y.data
        bd = dem_resample(de, grid_x, grid_y, var="adjusted", wet=True, flag=0, function="gauss")
        #        m['depth'].loc[dict(nSCHISM_hgrid_node=indices_of_nodes_in_bbox[ids])] =  - bd

        m["depth"].loc[dict(nSCHISM_hgrid_node=indices_of_nodes_in_bbox)] = -bd


def fill_dem_on_mesh(mesh_file, dem_files, coastlines):

    m0 = pm.set(type="schism", mesh_file=mesh_file)
    m_ = m0.Dataset

    for dfile in dem_files:

        dem_ = pdem.Dem(dem_source=dfile)
        dem = dem_.Dataset

        ilats = dem.elevation.chunk("auto").chunks[0]
        ilons = dem.elevation.chunk("auto").chunks[1]

        idx = [sum(ilons[:i]) for i in range(len(ilons) + 1)]
        jdx = [sum(ilats[:i]) for i in range(len(ilats) + 1)]

        blon = list(zip(idx[:-1], idx[1:]))
        blat = list(zip(jdx[:-1], jdx[1:]))

        perms = [(x, y) for x in blon for y in blat]

        populate(dem, perms, m_, ne_i, buffer=5)

        print(f"done with {dfile}")
