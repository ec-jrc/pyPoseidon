Framework for Sea Level Hydrodynamic simulations
================================================

[![Documentation Status](https://readthedocs.org/projects/pyposeidon/badge/?version=latest)](https://pyposeidon.readthedocs.io/en/latest/?badge=latest) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/ec-jrc/pyPoseidon) ![CI](https://github.com/ec-jrc/pyPoseidon/actions/workflows/conda_pip.yml/badge.svg) ![CI](https://github.com/ec-jrc/pyPoseidon/actions/workflows/conda_only.yml/badge.svg) ![CI](https://github.com/ec-jrc/pyPoseidon/actions/workflows/code_quality.yml/badge.svg) [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/ec-jrc/pyPoseidon/master?urlpath=%2Flab)

This is a development project utilising multiple solvers (currently `DELFT3D` & `SCHISM`) for simulating sea level height (currently only storm surge). The purpose is to create a simple, portable and transparent way of setting up, running and analysing hydrodynamic computations through python scripts and Jupyter Notebooks (http://jupyter.org). See Notebooks in Tutorial/ for relevant prototypes.

## Installation

There are two packages available. The base package which can be installed with

`conda install -c conda-forge pyposeidon-base`

while the full version which includes also the visualization dependencies needs

`conda install -c conda-forge pyposeidon`


### Prerequisities

`DELFT3D` needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. See Wiki for more details.

`SCHISM` needs to be compiled for your system. You can download it from  http://columbia.vims.edu/schism/tags/. See http://ccrm.vims.edu/schismweb/ for more info.


You can also install the solvers easily with conda

`conda install -c gbrey pschism delft3d4`


## Tests

There are several sets of tests. You can run `pyposeidon` unitests with

`pytest`

In order to test also the solver integration use

`pytest --runschism`

or

`pytest --rundelft`

if you are using a local installation of the solvers please specify the `PATH` to the executables in your system such as

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
