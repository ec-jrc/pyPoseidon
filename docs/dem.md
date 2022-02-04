<style>body {text-align: justify}</style>

The dem module handles the pre-processing of bathymetric/topographic data.


### Setup

- extent

The geometry's extent could be provided. In the most simple case that is a `lat/lon` box that defines the area of interest. Otherwise the full dataset will be used. This might be problematic for large files e.g. global datasets, especially when the file is retrieved from a url (see below).

Without loss of generality we select below *Iceland* as a test case.

```python
# define in a dictionary the properties of the model..
extend = {
    "lon_min": -25.0,  # lat/lon window
    "lon_max": -9.0,
    "lat_min": 56.0,
    "lat_max": 74.0,
}
```

- source

The `Dem` class supports all the functionality of `xarray`. Thus the `dem_source` argument can be:

#### Local file (or list of files) :

This could easily be defined as e.g.

```python
source = '/path/to/file/dem.nc'
```

A number of formats are supported, namely `geotiff, grib, netcdf, zarr`, etc.

#### URL 

For relative small areas an `erddap` server could be used e.g.

```python
url = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm30plus"
```

!!! note

	This option would work only for small extents. For large areas (continental/global) using a previously locally stored file is advised. 

- Coastlines

The resolution of dem datasets is usually lower that the corresponding ones for coastlines.

Since the coastlines are the boundary for the mesh generation, it makes sense to have them matched to the bathymetric data.

`pyposeidon` tries to facilitate that (by default) when you provide a coastline dataset as argument. This can be negated by setting `adjust_dem = False` in the arguments.

There are a number of available datasets for coastlines. One option is to use the `cartopy` features as:

```python
# use cartopy to get coastlines
import cartopy.feature as cf
import geopandas as gp

cr = "i"

coast = cf.NaturalEarthFeature(
    category="physical",
    name="land",
    scale="{}m".format({"l": 110, "i": 50, "h": 10}[cr]),
)

ne_i = gp.GeoDataFrame(geometry=[x for x in coast.geometries()])
```
 
for [natural earth](https://www.naturalearthdata.com) or 

 ```python
 # use cartopy to get coastlines
 gi = cf.GSHHSFeature(scale="intermediate", levels=[1])
 iGSHHS = gp.GeoDataFrame(geometry=[x for x in gi.geometries()])
 ```

for the [GSHHG](http://www.soest.hawaii.edu/pwessel/gshhg/) dataset. 

The coastlines argument could also be a `shapefile`. Generally, whatever format `geopandas` reads should work.


### Retrieve dem Dataset

We can combine the info into a dictionary as

```python
dic = {"geometry":extend, "dem_source":source, "coastlines":iGSHHS}
```

and retrieve the Dataset with

```python
import pyposeidon.dem as pdem
d = pdem.Dem(**dic)
```


!!! note

	`xarray` is using `dask` for a lazy read and all data will be loaded into memory when needed. 


### Resample on Mesh

In order to create a model, bathymetric data need to be interpolated on the corresponding mesh. If such a mesh is given, bathymetric data on the `x,y` nodes can be provided as 

```python
d = pdem.Dem(**dic, grid_x = x, grid_y = y)
```

 
### Output to file
 
 
Once the dataset is produced, it can be stored in the solver's appropriate format. This is relevant for solvers that read the bathymetry in a separate file (as is the case for `D3D`).
 
```python
pdem.to_output(d.Dataset,solver_name='d3d', rpath='./test/')
```

