<style>body {text-align: justify}</style>

`pyposeidon` is designed to be modular so that users can use any available functionality independently.

These modules include:

### Model

A top level function (see [model](api.md#pyposeidon.model)) that defines the model's characteristics, logging and most notably provides a `read_model()` option for reproducing an existing model by parsing a `json` file that defines the required settings.

### Boundary

This module handles mesh boundary extraction from coastlines or user provided shapefiles. It underpins the *Mesh* module.

### D3d

The wrapper of the `delft3d` solver is defined within a single file (see [d3d]())

### Schism 

The wrapper of the `schism` solver (see [schism](api.md#pyposeidon.schism)).

### Mesh

This is the module that handles the generation or incorporation of a given mesh (see [mesh](api.md#pyposeidon.mesh). It includes support for two mesh generators, [Jigsaw](https://github.com/dengwirda/jigsaw) and [GMSH](http://gmsh.info).

### Mjigsaw

Although `jigsaw` now provides a python API this module precedes that and instead wraps around the C++ binary.  


### Mgmsh

This handles the integration of both the binary and python bindings of the `gmsh` mesh generator. 

### Dem

This is the module that handles the pre-processing of bathymetric data (see [Dem](api.md#pyposeidon.dem)). It utilizes `xarray` to read and `pyresample` to interpolate on mesh points. In addition, has the ability to match coastlines to bathymetry.

### Meteo

A number of atmospheric data could be required for a geoflow simulation. These data could be in a variety of formats (`netcdf, grib, geotiff`, etc.). The meteo module uses `xarray` to read them into a dataset. This includes combining multiple files into a single dataset in terms of the time reference variable. 

### Tide

This module handles the tidal configuration. *Under development*.

### Utils

Finally, there are a number of supporting functions that are used by the main modules above that could also provide additional functionality on a stand-alone way.

- *pplot* : Provides an accessor to `xarray` for visualisations of Schism's unstructured mesh based on `matplolib`

- *hplot* : Provides an accessor to `xarray` for visualisations of Schism's unstructured mesh based on `holoviews`.

- *mplot* : Provides an accessor to `xarray` for visualisations of Schism's unstructured mesh based on `mayavi`.

- *seam* : Modify a global mesh so that it can be viewed in 2D without distortions.

- *cast* : Supports the hindcast or re-forecasting procedure.
