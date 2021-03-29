## Structure

pyPoseidon is designed to be modular so that users can use any available functionality independently.

These modules are:

### model

A level function (see [model]()) that defines the model characteristics, logging and most notably provides a read_model() option for reading in an existing model in the form of .json file that defines the required attributes.

### geometry

### d3d

The wrapper of the delft3d solver is defined within a single file (see [d3d]())

### schism 

The wrapper of the schism solver (see [schism](50_API.md#schism)).

### grid

This is the module that handles the generation or incorporation of a given mesh (see [grid.py]()). It includes support for 2 grid generators, [Jigsaw](https://github.com/dengwirda/jigsaw) and [GMSH](http://gmsh.info).

### jigsaw

Although Jigsaw now provides a python API this module precedes that and instead wraps around the C++ jigsaw solver.  


### ugmsh

GMSH has also a python API and this is used here for integration with pyPoseidon.

### dem

This is the module that handles the pre-processing of bathymetric data (see [dem]()). It utilizes xarray to read and pyresample to resample on mesh points. In addition, has the ability to match coastlines to bathymetry, again using pyresample and appropriate masking. See [adjust]() for more info.


### meteo

A number of atmospheric data could be required for a geoflow simulation. These data could be in a variety of formats (netcdf, grib, geotiff, etc.). The meteo module uses xarray to read them into a dataset. This includes combining multiple files into a single dataset in terms of the time reference variable. 

### tide

This module handles the solver's tidal configuration. *Under development*.

### utils

Finally, there are a number of supporting functions that are used by the main modules above that could also provide additional functionality on a stand-alone way.

- *pplot* : Provides an accessor to xarray for visualizations of schism's unstructured mesh based on matplolib

- *hplot* : Provides an accessor to xarray for visualizations of schism's unstructured mesh based on holoviews.

- *vplot* : Provides an accessor to xarray for visualizations of schism's unstructured mesh based on mayavi.

- *fix* : Adjusts bathymetric data to a given coastline.

- *cast* : Supports the hindcast or re-forecasting procedure.

- *seam* : Modify a global mesh so that it can be viewed in 2D without distortions.

