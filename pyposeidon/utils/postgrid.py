#!/usr/bin/env python
# coding: utf-8

import os

import numpy as np
import pandas as pd
import geopandas as gp
import pygeos
import xarray as xr
from tqdm import tqdm
import shapely


def nreduce(dnodes, delems, grid, elems, tag="tag"):

    ### Drop nodes

    grid = grid.drop(dnodes)  # drop nodes
    grid = grid.rename_axis(tag).reset_index()

    ### drop elements
    dropels = np.unique(delems)

    elems = elems.drop(dropels)
    elems.reset_index(inplace=True, drop=True)

    ### Re-index

    A, idxA = np.unique(grid[tag], return_inverse=True)
    B, idxB = np.unique(elems["a"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["a"] = idxA[IDX][idxB]
    B, idxB = np.unique(elems["b"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["b"] = idxA[IDX][idxB]
    B, idxB = np.unique(elems["c"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["c"] = idxA[IDX][idxB]

    return grid, elems


def drop(nodes, elems, bnodes, dq):

    nodes = nodes.drop(dq)  # drop nodes
    nodes = nodes.rename_axis("tag").reset_index()  # reset index

    ### Re-index tessalation

    A, idxA = np.unique(nodes["tag"], return_inverse=True)
    B, idxB = np.unique(elems["a"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["a"] = idxA[IDX][idxB]
    B, idxB = np.unique(elems["b"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["b"] = idxA[IDX][idxB]
    B, idxB = np.unique(elems["c"], return_inverse=True)
    IDX = np.in1d(A, B)
    elems["c"] = idxA[IDX][idxB]

    bnodes = bnodes.drop(bnodes.loc[bnodes.node.isin(dq)].index)
    sg = nodes.reset_index().set_index("tag")
    idx = bnodes.node.values, "index"
    nvs = sg.loc[idx].values
    bnodes.node = nvs

    bnodes = bnodes.reset_index(drop=True)
    bnodes.index.name = "bnodes"

    nodes = nodes.drop("tag", axis=1)

    return nodes, elems, bnodes


def check(g, shp, bad):

    # ## Read Grid
    #    g = pg.grid(type='tri2d',grid_file='/eos/jeodpp/data/projects/FLOODS-COAST/floods-coast/OPER/grids/eur_fixed.gr3')

    lon_min = g.Dataset.SCHISM_hgrid_node_x.values.min()
    lon_max = g.Dataset.SCHISM_hgrid_node_x.values.max()
    lat_min = g.Dataset.SCHISM_hgrid_node_y.values.min()
    lat_max = g.Dataset.SCHISM_hgrid_node_y.values.max()

    c = shp.cx[lon_min:lon_max, lat_min:lat_max]

    # ## Test polygons

    d = g.Dataset

    x = d.SCHISM_hgrid_node_x.values
    y = d.SCHISM_hgrid_node_y.values
    tri = d.SCHISM_hgrid_face_nodes.values

    nodes = pd.DataFrame({"lon": x, "lat": y})

    elems = pd.DataFrame(tri, columns=["a", "b", "c"])

    bnodes = g.Dataset[["node", "id", "type"]].to_dataframe()

    # drop bad
    dpoints = [bad]
    # get corresponding elements
    ma = elems.a.isin(dpoints)
    mb = elems.b.isin(dpoints)
    mc = elems.c.isin(dpoints)
    ids1 = elems[ma].index
    ids2 = elems[mb].index
    ids3 = elems[mc].index

    dropels = list(
        np.hstack([ids1, ids2, ids3])
    )  # these are the elements to be dropped - the ones which have dropped nodes

    dropels = np.unique(dropels)

    # take care of the boundary nodes
    dns = np.unique(elems.loc[dropels, ["a", "b", "c"]].values.flatten())
    # get boundary nodes
    ibs = [x for x in dns if x not in dpoints]

    maxb = bnodes.id.min()
    # check if new boundary nodes merge
    if bnodes.node.isin(ibs).sum() > 0:
        if bnodes.node.isin(ibs).sum() > 0:
            ids_ = np.unique(bnodes.loc[bnodes.node.isin(ibs), "id"].values)
            ids_.sort()
            if ids_.shape[0] > 1:
                itype = bnodes.loc[bnodes.id == ids_[-1], ["type"]].values[0]
                bnodes.loc[bnodes.id.isin(ids_[:-1]), "id"] = ids_[-1]
                ibs = pd.DataFrame(
                    {"node": ibs, "type": itype, "id": ids_[-1]},
                    index=np.arange(len(ibs)),
                )
            else:
                itype, iid = bnodes.loc[bnodes.node.isin(ibs), ["type", "id"]].values[0]
                ibs = pd.DataFrame({"node": ibs, "type": itype, "id": iid}, index=np.arange(len(ibs)))
    else:
        maxb -= 1
        ibs = pd.DataFrame({"node": ibs, "type": 1, "id": maxb}, index=np.arange(len(ibs)))

    bnodes = pd.concat([bnodes, ibs], ignore_index=True)

    bnodes = bnodes.drop_duplicates()

    # remove
    nodes, elems = nreduce(dpoints, dropels, nodes, elems)

    bnodes = bnodes.drop(bnodes.loc[bnodes.node.isin(dpoints)].index)
    sg = nodes.reset_index().set_index("tag")
    idx = bnodes.node.values, "index"
    nvs = sg.loc[idx].values
    bnodes.node = nvs

    bnodes = bnodes.reset_index(drop=True)

    bnodes.index.name = "bnodes"

    bnodes = bnodes.drop_duplicates("node")

    nodes = nodes.drop("tag", axis=1)

    # #### Check for hanging nodes

    # Look for hanging nodes
    tri3 = elems.values[:, :3]
    q = np.unique(tri3.flatten())  # all the unique nodes in elements

    dq = list(set(range(nodes.shape[0])) - set(q))  # the ones that are in gcrop but not in elems

    dq.sort()

    if len(dq) > 0:
        nodes, elems, bnodes = drop(nodes, elems, bnodes, dq)

    # ### Find the invalid nodes (that cross the coasts)
    cos = pygeos.from_shapely(c.geometry)
    cos_ = pygeos.set_operations.union_all(cos)

    gps = pygeos.points(list(nodes.values))

    gtree = pygeos.STRtree(gps)

    invs = gtree.query(cos_, predicate="contains").tolist()

    # ### Find invalid elements (that cross land)

    # cells to polygons
    ap = nodes.loc[elems.a]
    bp = nodes.loc[elems.b]
    cp = nodes.loc[elems.c]

    elems["ap"] = ap.values.tolist()
    elems["bp"] = bp.values.tolist()
    elems["cp"] = cp.values.tolist()

    n = 2
    al = elems.ap + elems.bp + elems.cp + elems.ap
    coords = [[l[i : i + n] for i in range(0, len(l), n)] for l in al]
    elems["coordinates"] = coords

    jig = pygeos.polygons(coords)

    jtree = pygeos.STRtree(jig)

    jig_ = pygeos.set_operations.union_all(jig)

    cross = pygeos.set_operations.intersection(jig_, cos_)

    # #### convert to dataframe

    fd = pd.DataFrame({"overlap": pygeos.to_wkt(cross)}, index=[0])

    fd["overlap"] = fd["overlap"].apply(shapely.wkt.loads)

    gover = gp.GeoDataFrame(fd, geometry="overlap")

    # #### Reject small injuctions
    ipols = gover.explode(index_parts=True).loc[0]

    ipols.columns = ["geometry"]

    mask = ipols.area.values == 0.0

    ipols = ipols[~mask].reset_index(drop=True)
    ipols = gp.GeoDataFrame(ipols)

    # Sometimes the edges cross the coastlines but the nodes not

    # ## Get overlapping elements for each contour

    jcos = pygeos.from_shapely(ipols.geometry)

    maxb = bnodes.id.min()

    if len(jcos) > 0:
        m = 0
        for l in tqdm(range(len(jcos))):
            con = jcos[l]
            m += 1
            jig = pygeos.polygons(elems.coordinates.to_list())
            jtree = pygeos.STRtree(jig)
            ci = jtree.query(con, predicate="intersects").tolist()
            if not ci:
                continue
            delems = elems.loc[ci].index.values
            dnodes = np.unique(elems.loc[ci, ["a", "b", "c"]].values.flatten())
            #    print (ci, dnodes, delems)

            inp = pygeos.points(list(nodes.loc[dnodes].values))
            itree = pygeos.STRtree(inp)
            ipoints = itree.query(cos_, predicate="contains").tolist()
            ipoints = nodes.loc[dnodes].index[ipoints].values

            # find elements that use these points
            ma = elems.a.isin(ipoints)
            mb = elems.b.isin(ipoints)
            mc = elems.c.isin(ipoints)
            ids1 = elems[ma].index
            ids2 = elems[mb].index
            ids3 = elems[mc].index

            dropels = list(
                np.hstack([ids1, ids2, ids3])
            )  # these are the elements to be dropped - the ones which have dropped nodes

            dropels = np.unique(dropels)

            dds = np.unique(np.concatenate((delems, dropels), axis=0))

            # get all nodes
            dns = np.unique(elems.loc[dds, ["a", "b", "c"]].values.flatten())
            # get boundary nodes
            ibs = [x for x in dns if x not in ipoints]

            #    print (ipoints, ibs)

            # check if new boundary nodes merge
            if bnodes.node.isin(ibs).sum() > 0:
                ids_ = np.unique(bnodes.loc[bnodes.node.isin(ibs), "id"].values)
                ids_.sort()
                if ids_.shape[0] > 1:
                    itype = bnodes.loc[bnodes.id == ids_[-1], ["type"]].values[0]
                    bnodes.loc[bnodes.id.isin(ids_[:-1]), "id"] = ids_[-1]
                    ibs = pd.DataFrame(
                        {"node": ibs, "type": itype, "id": ids_[-1]},
                        index=np.arange(len(ibs)),
                    )
                else:
                    itype, iid = bnodes.loc[bnodes.node.isin(ibs), ["type", "id"]].values[0]
                    ibs = pd.DataFrame(
                        {"node": ibs, "type": itype, "id": iid},
                        index=np.arange(len(ibs)),
                    )
            else:
                maxb -= 1
                ibs = pd.DataFrame({"node": ibs, "type": 1, "id": maxb}, index=np.arange(len(ibs)))

            bnodes = pd.concat([bnodes, ibs], ignore_index=True)

            bnodes = bnodes.drop_duplicates()

            nodes, elems = nreduce(ipoints, dds, nodes, elems)

            # Renumber boundary nodes
            if ipoints.size > 0:
                bnodes = bnodes.drop(bnodes.loc[bnodes.node.isin(ipoints)].index)
                sg = nodes.reset_index().set_index("tag")
                #        print(ipoints.size)
                idx = bnodes.node.values, "index"
                nvs = sg.loc[idx].values
                bnodes.node = nvs

            nodes = nodes.drop("tag", axis=1)

    # #### Check for hanging nodes

    # Look for hanging nodes
    tri3 = elems.values[:, :3]
    q = np.unique(tri3.flatten())  # all the unique nodes in elements

    dq = list(set(range(nodes.shape[0])) - set(q))  # the ones that are in gcrop but not in elems

    dq.sort()

    if len(dq) > 0:
        nodes, elems, bnodes = drop(nodes, elems, bnodes, dq)

    # Get the largest continous area
    jels = pygeos.polygons(elems.coordinates.values.tolist())
    wat = pygeos.set_operations.coverage_union_all(jels)
    w = pd.DataFrame({"overlap": pygeos.to_wkt(wat)}, index=[0])
    w["overlap"] = w["overlap"].apply(shapely.wkt.loads)
    gw = gp.GeoDataFrame(w, geometry="overlap")

    gw = gw.explode(index_parts=True)
    gw = gw.droplevel(0)
    gw.columns = ["geometry"]
    gw["length"] = gw["geometry"][:].length
    gw = gw.sort_values(by="length", ascending=0)  # optional
    gw = gw.reset_index(drop=True)

    # indentify the elements of the large polygon
    cgw = pygeos.from_shapely(gw.loc[1:].geometry)
    cgw_ = pygeos.set_operations.union_all(cgw)
    jtree_ = pygeos.STRtree(jels)
    invs = jtree_.query(cgw_, predicate="intersects").tolist()

    if len(invs) > 0:
        # Sort the elements (some shouldn't be there)

        qnodes = np.unique([elems.loc[x, ["a", "b", "c"]].values.astype(int) for x in invs])
        nnodes = pygeos.points(list(nodes.loc[qnodes].values))
        ntree = pygeos.STRtree(nnodes)
        nels_ = pygeos.from_shapely(gw.loc[0].geometry.buffer(0.00001))
        nevs = ntree.query(nels_, predicate="intersects").tolist()
        pns = qnodes[nevs]
        dpoints = [x for x in qnodes if x not in pns]

        # find elements that use these points
        ma = elems.a.isin(dpoints)
        mb = elems.b.isin(dpoints)
        mc = elems.c.isin(dpoints)
        ids1 = elems[ma].index
        ids2 = elems[mb].index
        ids3 = elems[mc].index

        dropels = list(
            np.hstack([ids1, ids2, ids3])
        )  # these are the elements to be dropped - the ones which have dropped nodes

        dropels = np.unique(dropels)

        # Remove the invalid elements
        nodes, elems = nreduce(dpoints, dropels, nodes, elems)

        # reindex boundary nodes
        bnodes = bnodes.drop(bnodes.loc[bnodes.node.isin(dpoints)].index)
        sg = nodes.reset_index().set_index("tag")
        idx = bnodes.node.values, "index"
        nvs = sg.loc[idx].values
        bnodes.node = nvs

        bnodes = bnodes.reset_index(drop=True)

        bnodes.index.name = "bnodes"

        nodes = nodes.drop("tag", axis=1)

    # Look for hanging nodes
    tri3 = elems.values[:, :3]
    q = np.unique(tri3.flatten())  # all the unique nodes in elements

    dq = list(set(range(nodes.shape[0])) - set(q))  # the ones that are in gcrop but not in elems

    dq.sort()

    if len(dq) > 0:
        nodes, elems, bnodes = drop(nodes, elems, bnodes, dq)

    # check if empty boundary

    lk = -1
    for k in range(-1, bnodes.id.min().astype(int) - 1, -1):
        if not bnodes.loc[bnodes.id == k].empty:
            bnodes.loc[bnodes.id == k, "id"] = lk
            lk -= 1

    bnodes = bnodes.drop_duplicates("node")

    # ### create the new dataset

    nod = (
        nodes.loc[:, ["lon", "lat"]]
        .to_xarray()
        .rename(
            {
                "index": "nSCHISM_hgrid_node",
                "lon": "SCHISM_hgrid_node_x",
                "lat": "SCHISM_hgrid_node_y",
            }
        )
    )
    nod = nod.drop_vars("nSCHISM_hgrid_node")

    depx = xr.Dataset({"depth": (["nSCHISM_hgrid_node"], np.zeros(nod.nSCHISM_hgrid_node.shape[0]))})

    elsx = xr.DataArray(
        elems.loc[:, ["a", "b", "c"]].values,
        dims=["nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"],
        name="SCHISM_hgrid_face_nodes",
    )

    ngr = xr.merge([nod, depx, elsx, bnodes.to_xarray()])  # total

    return ngr
