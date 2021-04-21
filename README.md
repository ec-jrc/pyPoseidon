Framework for Sea Level Hydrodynamic simulations
================================================

[![Documentation Status](https://readthedocs.org/projects/pyposeidon/badge/?version=latest)](https://pyposeidon.readthedocs.io/en/latest/?badge=latest) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/brey/pyPoseidon) ![CI](https://github.com/brey/pyPoseidon/actions/workflows/conda_and_nested_venv.yml/badge.svg) ![CI](https://github.com/brey/pyPoseidon/actions/workflows/conda.yml/badge.svg) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/brey/pyPoseidon/master?urlpath=%2Flab)

This is a development project utilising multiple solvers (currently DELFT3D & SCHISM) for simulating sea level height (currently only storm surge). The purpose is to create a simple, portable and transparent way of setting up, running and analysing hydrodynamic computations through python scripts and Jupyter Notebooks (http://jupyter.org). See Notebooks in Tutorial/ for relevant prototypes.

## Installation


`conda install -c gbrey pyposeidon`

Afterwards, for now, one needs to install gmsh manually with

`pip install gmsh`

**Note**: Due to an upstream issue, *pydap* needs to be installed manually. See *environment.yml* for info.

### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. See Wiki for more details.

SCHISM needs to be compiled for your system. You can download it from  http://columbia.vims.edu/schism/tags/. See http://ccrm.vims.edu/schismweb/ for more info.


You can also install the solvers easily with conda

`conda install -c gbrey pschism delft3d4`


## Tests

There are several sets of tests. You can run pyPoseidon unitests with

`pytest`

In order to test also the solver integration use

`pytest --runschism`

or

`python --rundelft`

if you are using a local installation of the solvers please specify the PATH to the executables in your system such as

`export D3D = '/path_to_folder_bin/lnx64/flow2d3d/'`

`export LD3D = '/path_to_folder_bin/lnx64/flow2d3d/'`

`export SCHISM = '/path_to_schism_executable'`

## docs

```
mkdocs build
mkdocs serve
```

## License
* The project is released under the EUPL v1.2 license.
