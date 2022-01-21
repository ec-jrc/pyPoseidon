<style>body {text-align: justify}</style>

There are a number of ways to install `pyposeidon`.

### conda

`conda install -c conda-forge pyposeidon`

**Note**: Due to an upstream issue, `pydap` needs to be installed manually. See `dependencies/python_deps10.yml` for info.

### pip

For using `pip` please note that your system needs to have:

- `python>=3.8`
- `geos`
- `gdal=3.2.1`
- `proj<8`

Then 

`pip install pyposeidon`


### virtualenv

You could also use `conda` for the above system dependencies by creating a dedicated env with 

`conda create -n pyPoseidon pip python=3.8 geos gdal=3.2.1 proj=7 poetry`

Afterwards, activate the new `conda` environment, create a `virtualenv` and install the dependencies using `poetry`:

`conda activate pyPoseidon`

`python3 -m venv .venv`

`source .venv/bin/activate`

`poetry install`


You are ready to go!


### Solver integration

DELFT3D needs to be compiled for your system. You can download it from [deltares](http://oss.deltares.nl/web/delft3d/source-code).

SCHISM needs to be compiled for your system. You can download it from  [github](https://github.com/schism-dev/schism).


Alternatively, you can also install the solvers with conda

`conda install -c gbrey pschism delft3d4`


### Tests

You can use `pytest` to run the tests. Since the required datasets are substantial in size these data will be downloaded once you initiate `pytest` for the first time. 

So, install `pytest` either via 

`conda install pytest` 

if you are using a `conda` based env or by 


`pip install pytest`

if you are using `venv`.

Then


`pytest`

In order to test also the solver integration use

`pytest --runschism`

or

`pytest --rundelft`

if you are using a local installation of the solvers please specify the `PATH` to the executables in your system such as

`export D3D = '/path_to_folder_bin/lnx64/flow2d3d/'`

`export LD3D = '/path_to_folder_bin/lnx64/flow2d3d/'`

`export SCHISM = '/path_to_schism_executable'`

