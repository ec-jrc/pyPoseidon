<style>body {text-align: justify}</style>

The top module that handles the main setup of the instance.


### Setup

A number of parameters is provided by the user :

- *solver_name*:
	Select the solver to be used e.g. *schism*.
- *geometry*:
	In the most simple case that is a lat/lon box that defines the area of interest.
- *coastlines*:
	A coastlines *GeoDataFrame* or *shapefile* providing boundaries to the mesh.
- *mesh_generator*:
	Set the backend for creating a mesh e.g. `jigsaw`.
- *start_date*:
	timestamp that the simulation starts.
- *time_frame*:
	Duration of the simulation after the *start_date*. Alternatively one can set *end_date*.
- *meteo_source*:
	Source for atmospheric data for forcing the model.
- *dem_source*:
	Topobathy Dataset.

All necessary info can be incorporated into a dictionary as

```python
dic = {
    "solver_name": "schism",
    "geometry": {"lon_min": -25, "lon_max": -12.0, "lat_min": 56.0, "lat_max": 74.0},
    "coastlines": "/path/to/coastal/shapefile.shp",
    "mesh_generator": "jigsaw",
    "start_date": "2017-10-1 0:0:0",
    "time_frame": "12h",
    "meteo_source": ["/path/to/meteo/file.grib"],
    "dem_source": "./path/to/dem/file.nc",
}
```

Having all the attributes defined, the model can be initiated as

```python
import pyposeidon.model as pm
b = pm.set(**dic)
```

Executing the model can be done incrementally or in one step. The steps involved are:

```py
b.create() # constructs all required parts e.g. mesh, dem, meteo, etc.
b.output() # save to files
b.save() # saves the json model reference file
b.set_obs() # setup station points
b.run() # execute
```

or in one step as

```py
b.execute()
```

The various datasets incorporated into the model can be accessed independently as attributes, namely :

```py
b.meteo.Dataset # forcing
b.mesh.Dataset # mesh
b.dem.Dataset # bathymetry
b.coastlines # cooastlines used
```

### Output

The output of the simulations could be in separate files (due to `mpi`). These files can be integrated into a single Dataset with

```py
b.get_output_data() # integrate output
```

and the data are available as

```py
b.data.Dataset # output Dataset
```


### Read

If a model is created by `pyposeidon` there is a `json` file that describes the model. This file serves as a signature for the specific model. It can be used to easily reproduce the results by reading the file into `pyposeidon` and executing it as

```py
a = pyposeidon.model.read('./path/to/schism_model.json')
a.execute() #execute it again
```

When a model is created by other means, it can still be incorporated into `pyposeidon` with

```py
c = pm.set(solver_name='schism', rfolder = './path/to/folder/',
 	 		load_mesh=True, load_meteo=True)
c.mesh.Dataset
c.meteo.Dataset
```
