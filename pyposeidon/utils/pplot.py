"""
Visualization module based on matplotlib

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import geopandas as gp
import pandas as pd
import shapely
from pyposeidon.utils.quads2tr import quads_to_tris

import sys
import os

ffmpeg = sys.exec_prefix + "/bin/ffmpeg"
os.environ["FFMPEG_BINARY"] = ffmpeg
from matplotlib import animation


matplotlib.rc("animation", html="html5")
plt.rcParams["animation.html"] = "jshtml"
plt.rcParams["animation.embed_limit"] = "200."


def __init__(dark_background=False):
    # set plt style
    if dark_background:
        plt.style.use("dark_background")


@xr.register_dataset_accessor("gplot")
class gplot(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def contourf(self, ax=None, x=None, y=None, z=None, tname="time", **kwargs):
        if not ax:
            fig = plt.figure(figsize=(12, 8))

        if len(self._obj[x].shape) > 2:
            grid_x = self._obj[x].values[0, :, :]
            grid_y = self._obj[y].values[0, :, :]
        else:
            grid_x = self._obj[x].values
            grid_y = self._obj[y].values

        z_ = self._obj[z].values
        t = self._obj[tname].values

        vmin = kwargs.get("vmin", z_.min())
        vmax = kwargs.get("vmax", z_.max())

        nv = kwargs.get("nv", 10)

        title = kwargs.get("title", None)

        coasts = kwargs.get("coastlines", False)

        vrange = np.linspace(vmin, vmax, nv, endpoint=True)
        ## CHOOSE YOUR PROJECTION
        #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        ax = plt.axes(projection=ccrs.PlateCarree())
        # Limit the extent of the map to a small longitude/latitude range.

        ax.set_aspect("equal")
        ims = []
        for i in range(len(t)):
            im = ax.contourf(
                grid_x,
                grid_y,
                z_[i, :, :],
                vrange,
                vmin=vmin,
                vmax=vmax,
                transform=ccrs.PlateCarree(),
            )
            #        im = ax.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
            add_arts = im.collections
            text = "time={}".format(t[i])
            # te = ax.text(90, 90, text)
            an = ax.annotate(text, xy=(0.05, 1.05), xycoords="axes fraction")
            ims.append(add_arts + [an])
        if title:
            ax.set_title(title)
        # ax.set_global()
        if coasts:
            ax.coastlines("50m")
        ax.set_extent([grid_x.min(), grid_x.max(), grid_y.min(), grid_y.max()])
        ax.gridlines(draw_labels=True)

        # cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])
        cbar = fig.colorbar(im, ticks=vrange, orientation="vertical")  # ,fraction=0.046, pad=0.04)
        # plt.colorbar()

        v = animation.ArtistAnimation(fig, ims, interval=200, blit=False, repeat=False)

        plt.close()

        return v

    def update_quiver(num, Q, U, V, step):
        """updates the horizontal and vertical vector components by a
        fixed increment on each frame
        """

        Q.set_UVC(U[num, ::step, ::step], V[num, ::step, ::step])

        return (Q,)

    def quiver(self, ax=None, x=None, y=None, z=None, tname="time", **kwargs):
        U = self._obj[z].values[:, :, :, 0]
        V = self._obj[z].values[:, :, :, 1]

        if len(self._obj[x].shape) > 2:
            X = self._obj[x].values[0, :, :]
            Y = self._obj[y].values[0, :, :]
        else:
            X = self._obj[x].values
            Y = self._obj[y].values

        if not ax:
            fig = plt.figure(figsize=(12, 8))
            ax = plt.axes(projection=ccrs.PlateCarree())
        crs = ccrs.PlateCarree()
        ax.set_aspect("equal")

        land_50m = cfeature.NaturalEarthFeature(
            "physical",
            "land",
            "50m",
            edgecolor="face",
            facecolor=cfeature.COLORS["land"],
            zorder=0,
        )

        sea_50m = cfeature.NaturalEarthFeature(
            "physical",
            "ocean",
            "50m",
            edgecolor="face",
            facecolor=cfeature.COLORS["water"],
            zorder=0,
        )

        title = kwargs.get("title", None)

        ax.coastlines("50m")
        ax.add_feature(land_50m)
        ax.add_feature(sea_50m)

        scale = kwargs.get("scale", 1.0)  # change accordingly to fit your needs
        step = kwargs.get("step", 1)  # change accordingly to fit your needs

        Q = ax.quiver(
            X[::step, ::step],
            Y[::step, ::step],
            U[0, ::step, ::step],
            V[0, ::step, ::step],
            pivot="mid",
            color="k",
            angles="xy",
            scale_units="xy",
            scale=scale,
            transform=crs,
        )

        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())
        ax.set_title(title)
        # ax.set_global()

        plt.close()
        # you need to set blit=False, or the first set of arrows never gets
        # cleared on subsequent frames
        v = animation.FuncAnimation(
            fig,
            update_quiver,
            fargs=(Q, U, V, step),
            frames=range(0, np.size(t)),
            blit=False,
            repeat=False,
        )  # , interval=1)

        return v


def update_qframes(num, Q, U, V):
    """updates the horizontal and vertical vector components by a
    fixed increment on each frame
    """

    Q.set_UVC(U[num, :], V[num, :])

    return (Q,)


@xr.register_dataset_accessor("pplot")
# @xr.register_dataarray_accessor('pplot')


class pplot(object):
    def __init__(self, xarray_obj):
        self._obj = xarray_obj

    def contour(self, ax=None, it=None, **kwargs):
        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")
        t_var = kwargs.get("t", "time")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        try:
            t = self._obj[t_var][:].values
        except:
            pass

        tes = self._obj[tes_var].values[:, :4]

        # check tesselation

        if tes.max() == x.size:
            tes = tes - 1  # fortran/python conversion

        try:
            tes = tes.fillna(-1.0)  # adapt to schism >=5.10
            tes = tes.values
        except:
            pass  # given as numpy array

        var = kwargs.get("var", "depth")
        z = kwargs.get("z", self._obj[var].values[it, :].flatten())

        # sort out quads
        try:
            mask = tes[:, 3] == -1
            tr3 = tes[mask][:, :3]
            tr3_ = quads_to_tris(tes[~mask])
            if tr3_:
                tri3 = np.append(tr3, tr3_, axis=0).astype(int)
            else:
                tri3 = tr3.astype(int)
        except:
            tri3 = tes[:, :3].astype(int)

        if tri3.min() > 0:
            tri3 = tri3 - 1
            try:
                quads = quads - 1
            except:
                pass

        if not ax:
            fig, ax = plt.subplots(figsize=(12, 8))
        ## CHOOSE YOUR PROJECTION
        #   ax = plt.axes(projection=ccrs.Orthographic(x.mean(), y.mean()))
        #   ax = plt.axes(projection=ccrs.PlateCarree())
        #   ax.background_patch.set_facecolor('k')
        # ax = plt.axes()

        nv = kwargs.get("nv", 10)
        xy = kwargs.get("xy", (0.05, -0.1))
        title = kwargs.get("title", "contour plot for {}".format(var))

        # optional mask for the data
        mask = kwargs.get("mask", None)
        fv = kwargs.get("fv", -99999)
        if "mask" in kwargs:
            z = np.ma.masked_array(z, mask)
            z = z.filled(fill_value=fv)

        # check for nans
        imask = np.isnan(z)
        if imask.sum() > 0:
            z = np.ma.masked_array(z, imask)
            z = z.filled(fill_value=fv)

        vmin = kwargs.get("vmin", z.min())
        vmax = kwargs.get("vmax", z.max())
        vrange = np.linspace(vmin, vmax, nv, endpoint=True)

        for val in [
            "x",
            "y",
            "t",
            "it",
            "vmin",
            "vmax",
            "title",
            "nv",
            "tes",
            "mask",
            "xy",
            "z",
            "e",
            "fv",
            "var",
        ]:
            try:
                del kwargs[val]
            except:
                pass

        ax.set_aspect("equal")

        p = plt.tricontour(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, **kwargs)
        cbar = fig.colorbar(p, ticks=vrange, orientation="vertical")
        if it:
            text = "time={}".format(t[it])
            an = ax.annotate(text, xy=xy, xycoords="axes fraction")

        ax.set_title(title, pad=30)
        plt.xlabel("Longitude (degrees)")
        plt.ylabel("Latitude (degrees)")

        return p  # , ax

    def contourf(self, ax=None, it=None, **kwargs):
        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")
        t_var = kwargs.get("t", "time")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        try:
            t = self._obj[t_var][:].values
        except:
            pass

        tes = self._obj[tes_var].values[:, :4]

        # check tesselation

        if tes.max() == x.size:
            tes = tes - 1  # fortran/python conversion

        try:
            tes = tes.fillna(-1.0)  # adapt to schism >=5.10
            tes = tes.values
        except:
            pass  # given as numpy array

        # sort out quads
        try:
            mask = tes[:, 3] == -1
            tr3 = tes[mask][:, :3]
            tr3_ = quads_to_tris(tes[~mask])
            if tr3_:
                tri3 = np.append(tr3, tr3_, axis=0).astype(int)
            else:
                tri3 = tr3.astype(int)
        except:
            tri3 = tes[:, :3].astype(int)

        if tri3.min() > 0:
            tri3 = tri3 - 1
            try:
                quads = quads - 1
            except:
                pass

        var = kwargs.get("var", "depth")
        z = kwargs.get("z", self._obj[var].values[it, :].flatten())

        nv = kwargs.get("nv", 10)

        title = kwargs.get("title", "contourf plot for {}".format(var))

        ## CHOOSE YOUR PROJECTION
        #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        #    [fig,ax] = kwargs.get('figure',[plt.figure(figsize=(12,8)),plt.axes(projection=ccrs.PlateCarree())])
        #    ax.set_extent([x.min(), x.max(), y.min(), y.max()])
        #     ax.background_patch.set_facecolor('k')

        if not ax:
            fig = plt.figure(figsize=(12, 8))
            ax = plt.axes()

        # optional mask for the data
        mask = kwargs.get("mask", None)
        fv = kwargs.get("fv", -99999)
        if "mask" in kwargs:
            z = np.ma.masked_array(z, mask)
            z = z.filled(fill_value=fv)

        # check for nans
        imask = np.isnan(z)
        if imask.sum() > 0:
            z = np.ma.masked_array(z, imask)
            z = z.filled(fill_value=fv)

        vmin = kwargs.get("vmin", z.min())
        vmax = kwargs.get("vmax", z.max())

        vrange = np.linspace(vmin, vmax, nv, endpoint=True)

        xy = kwargs.get("xy", (0.05, -0.1))

        for val in [
            "x",
            "y",
            "t",
            "it",
            "z",
            "vmin",
            "vmax",
            "title",
            "nv",
            "tes",
            "mask",
            "xy",
            "var",
            "e",
            "fv",
            "figure",
        ]:
            try:
                del kwargs[val]
            except:
                pass

        ax.set_aspect("equal")

        p = ax.tricontourf(x, y, tri3, z, vrange, vmin=vmin, vmax=vmax, **kwargs)  # , transform=ccrs.PlateCarree() )
        cbar = plt.colorbar(p, ticks=vrange, orientation="vertical")
        if it:
            text = "time={}".format(t[it])
            an = ax.annotate(text, xy=xy, xycoords="axes fraction")

        ax.set_title(title, pad=30)
        ax.set_xlabel("Longitude (degrees)")
        ax.set_ylabel("Latitude (degrees)")

        return ax

    def quiver(self, ax=None, it=None, u=None, v=None, title=None, scale=0.1, color="k", **kwargs):
        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        t_var = kwargs.get("t", "time")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        try:
            t = self._obj[t_var][:].values
        except:
            pass

        if not ax:
            fig = plt.figure(figsize=(12, 8))
            ## CHOOSE YOUR PROJECTION
            #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
            #   ax = plt.axes(projection=ccrs.PlateCarree())
            #   ax.background_patch.set_facecolor('k')
            ax = plt.gca()

        title = kwargs.get("title", "vector plot for {}".format(title))
        xy = kwargs.get("xy", (0.05, -0.1))

        # optional mask for the data
        mask = kwargs.get("mask", None)
        fv = kwargs.get("fv", -99999)
        if "mask" in kwargs:
            u = np.ma.masked_array(u, mask)
            v = np.ma.masked_array(v, mask)
            v = v.filled(fill_value=fv)
            u = u.filled(fill_value=fv)

        for val in [
            "x",
            "y",
            "t",
            "it",
            "u",
            "v",
            "fv",
            "title",
            "tes",
            "xy",
            "scale",
            "mask",
            "color",
            "var",
        ]:
            try:
                del kwargs[val]
            except:
                pass

        ax.set_aspect("equal")

        p = plt.quiver(x, y, u, v, angles="xy", scale_units="xy", scale=scale, color=color, **kwargs)
        plt.xlabel("Longitude (degrees)")
        plt.ylabel("Latitude (degrees)")
        ax.set_title(title, pad=30)

        if it:
            text = "time={}".format(t[it])
            an = ax.annotate(text, xy=xy, xycoords="axes fraction")

        return ax

    def mesh(self, ax=None, lw: float = 0.5, markersize: float = 1.0, **kwargs):
        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        tes = self._obj[tes_var].values[:, :4]

        # check tesselation

        try:
            tes = tes.fillna(-1.0)  # adapt to schism >=5.10
            tes = tes.values
        except:
            pass  # given as numpy array

        # sort out quads
        try:
            mask = tes[:, 3] == -1
            tri3 = tes[mask][:, :3].astype(int)
            quads = tes[~mask].astype(int)
        except:
            tri3 = tes[:, :3].astype(int)
            quads = []

        if tri3.min() > 0:
            tri3 = tri3 - 1
            try:
                quads = quads - 1
            except:
                pass

        for val in ["x", "y", "e"]:
            try:
                del kwargs[val]
            except:
                pass

        if not ax:
            fig = plt.figure(figsize=(12, 8))
            ax = plt.gca()
            # ax = plt.axes(projection=ccrs.PlateCarree())
            # ax.background_patch.set_facecolor('k')

        # optionally plot boundary nodes
        bnodes = kwargs.get("bnodes", False)
        if bnodes:
            bs = []
            for group, ds in self._obj.groupby("id"):
                xi, yi = x[ds.node], y[ds.node]
                ib = shapely.geometry.MultiPoint(list(zip(xi, yi)))
                bs.append(ib)
            bls = gp.GeoDataFrame(geometry=bs)
            bls = bls[::-1].reset_index(drop=True)
            bls.plot(ax=ax, markersize=10, marker="x")
            del kwargs["bnodes"]

        ax.set_aspect("equal")
        g = plt.triplot(x, y, tri3, "go-", lw=lw, markersize=markersize, **kwargs)  # transform=ccrs.PlateCarree() )

        lw = kwargs.get("lw", plt.rcParams["lines.linewidth"])
        # https://stackoverflow.com/questions/52202014/how-can-i-plot-2d-fem-results-using-matplotlib
        for element in quads:
            x_ = [x[element[i]] for i in range(len(element))]
            y_ = [y[element[i]] for i in range(len(element))]
            plt.fill(x_, y_, edgecolor="green", fill=False, lw=lw)

        title = kwargs.get("title", "Mesh plot")
        ax.set_title(title, pad=30)
        ax.set_xlabel("Longitude (degrees)")
        ax.set_ylabel("Latitude (degrees)")

        return ax

    def qframes(self, ax=None, u=None, v=None, scale=0.01, color="k", **kwargs):
        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        t_var = kwargs.get("t", "time")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        t = self._obj[t_var][:].values

        cr = kwargs.get("coastlines", None)
        c_attrs = kwargs.get("coastlines_attrs", {})

        #        ax = plt.axes(projection=ccrs.PlateCarree())
        #  ax.set_extent([x.min(), x.max(), y.min(), y.max()])

        if not ax:
            fig = plt.figure(figsize=(12, 8))
            ax = plt.gca()

        ax.set_aspect("equal")

        title = kwargs.get("title", None)

        step = kwargs.get("step", 1)  # change accordingly to fit your needs

        Q = ax.quiver(
            x,
            y,
            u[0, :],
            v[0, :],
            pivot="mid",
            color=color,
            angles="xy",
            scale_units="xy",
            scale=scale,
        )

        #        if cr is not None:
        #            try:
        #                coastl = gp.GeoDataFrame.from_file(cr)
        #            except:
        #                coastl = gp.GeoDataFrame(cr)
        #            coastl.plot(ax=ax, **c_attrs)

        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_title(title)
        # ax.set_global()

        # you need to set blit=False, or the first set of arrows never gets
        # cleared on subsequent frames
        v = animation.FuncAnimation(
            fig,
            update_qframes,
            fargs=(Q, u, v),
            blit=False,
            repeat=False,
            frames=range(0, np.size(t)),
        )

        plt.close()

        return v

    def frames(self, **kwargs):
        cr = kwargs.get("coastlines", None)
        c_attrs = kwargs.get("coastlines_attrs", {})

        x_var = kwargs.get("x", "SCHISM_hgrid_node_x")
        y_var = kwargs.get("y", "SCHISM_hgrid_node_y")
        tes_var = kwargs.get("e", "SCHISM_hgrid_face_nodes")
        t_var = kwargs.get("t", "time")

        x = self._obj[x_var][:].values
        y = self._obj[y_var][:].values
        t = self._obj[t_var][:].values
        tes = self._obj[tes_var].values[:, :4]

        # check tesselation

        try:
            tes = tes.fillna(-1.0)  # adapt to schism >=5.10
            tes = tes.values
        except:
            pass  # given as numpy array

        # sort out quads
        try:
            mask = tes[:, 3] == -1
            tr3 = tes[mask][:, :3]
            tr3_ = quads_to_tris(tes[~mask])
            if tr3_:
                tri3 = np.append(tr3, tr3_, axis=0).astype(int)
            else:
                tri3 = tr3.astype(int)
        except:
            tri3 = tes[:, :3].astype(int)

        if tri3.min() > 0:
            tri3 = tri3 - 1  # fortran/python conversion

        var = kwargs.get("var", "depth")
        z = kwargs.get("z", self._obj[var].values)

        # set figure size
        xr = x.max() - x.min()
        yr = y.max() - y.min()
        ratio = yr / xr
        xf = 12
        yf = np.ceil(12 * ratio).astype(int)

        fig = plt.figure(figsize=(xf, yf))

        nv = kwargs.get("nv", 10)

        title = kwargs.get("title", None)

        # optional mask for the data
        mask = kwargs.get("mask", None)
        fv = kwargs.get("fv", -99999)
        if "mask" in kwargs:
            z = np.ma.masked_array(z, mask)
            z = z.filled(fill_value=fv)

        # check for nans
        imask = np.isnan(z)
        if imask.sum() > 0:
            z = np.ma.masked_array(z, imask)
            z = z.filled(fill_value=fv)

        vmin = kwargs.get("vmin", z.min())
        vmax = kwargs.get("vmax", z.max())
        vrange = np.linspace(vmin, vmax, nv, endpoint=True)

        ## CHOOSE YOUR PROJECTION
        #   ax = plt.axes(projection=ccrs.Orthographic(grid_x.mean(), grid_y.mean()))
        #  ax = plt.axes(projection=ccrs.PlateCarree())
        #   ax.background_patch.set_facecolor('k')
        # Limit the extent of the map to a small longitude/latitude range.
        #    ax.set_extent([x.min(), x.max(), y.min(), y.max()])

        ax = plt.axes()
        ax.set_aspect("equal")

        ims = []
        for i in range(len(t)):
            im = ax.tricontourf(x, y, tri3, z[i, :], vrange, vmin=vmin, vmax=vmax)  # , transform=ccrs.PlateCarree())
            #        im = ax.contourf(x,y,z[i,:,:],v,vmin=v1,vmax=v2,latlon=True)
            add_arts = im.collections
            text = "time={}".format(t[i])
            # te = ax.text(90, 90, text)
            an = ax.annotate(text, xy=(0.05, -0.1), xycoords="axes fraction")
            ims.append(add_arts + [an])

        #            if cr is not None: TO DO
        #                try:
        #                    coastl = gp.GeoDataFrame.from_file(cr)
        #                except:
        #                    coastl = gp.GeoDataFrame(cr)
        #                coastl.plot(ax=ax, **c_attrs)

        if title:
            ax.set_title(title)
        # ax.set_global()
        # ax.coastlines('50m')

        # cbar_ax = fig.add_axes([0.05, 0.05, 0.85, 0.05])
        cbar = fig.colorbar(im, ticks=vrange, orientation="vertical", fraction=0.017, pad=0.04)
        # plt.colorbar()

        v = animation.ArtistAnimation(fig, ims, interval=200, blit=False, repeat=False)

        plt.close()

        return v


# https://pandas.pydata.org/pandas-docs/stable/development/extending.html
@pd.api.extensions.register_dataframe_accessor("geo")
class GeoAccessor:
    def __init__(self, geopandas_obj):
        dic = {}
        tag = 0
        for l in range(geopandas_obj.shape[0]):
            lon = []
            lat = []
            try:
                contour = geopandas_obj.loc[l].geometry.boundary
            except:
                contour = geopandas_obj.loc[l].geometry

            if contour.type == "LineString":
                for x, y in contour.coords[:]:
                    lon.append(x)
                    lat.append(y)
                dic.update({"line{}".format(tag): {"longitude": lon, "latitude": lat}})
                tag += 1
            elif contour.type == "MultiLineString":
                for m in range(len(contour)):
                    lon = []
                    lat = []
                    for x, y in contour[l].coords[:]:
                        lon.append(x)
                        lat.append(y)
                    dic.update({"line{}".format(tag): {"longitude": lon, "latitude": lat}})
                    tag += 1

        dict_of_df = {k: pd.DataFrame(v) for k, v in dic.items()}
        df = pd.concat(dict_of_df, axis=0)
        df = df.drop_duplicates()  # drop the repeat value on closed boundaries
        self._validate(df)
        self._obj = df

    @staticmethod
    def _validate(obj):
        # verify there is a column latitude and a column longitude
        if "latitude" not in obj.columns or "longitude" not in obj.columns:
            raise AttributeError("Must have 'latitude' and 'longitude'.")

    @property
    def center(self):
        # return the geographic center point of this DataFrame
        lat = self._obj.latitude
        lon = self._obj.longitude
        return (float(lon.mean()), float(lat.mean()))

    def plot(self):
        # plot this array's data on a map, e.g., using Cartopy
        self._obj.plot.scatter(x="longitude", y="latitude")
        pass
