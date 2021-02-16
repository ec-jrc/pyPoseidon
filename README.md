Framework for Sea Level Hydrodynamic simulations
================================================

[![Build Status](https://dev.azure.com/breyiannis/pyPoseidon/_apis/build/status/brey.pyPoseidon?branchName=master)](https://dev.azure.com/breyiannis/pyPoseidon/_build/latest?definitionId=1&branchName=master) ![Azure DevOps coverage](https://img.shields.io/azure-devops/coverage/breyiannis/pyPoseidon/1) ![GitHub release (latest by date)](https://img.shields.io/github/v/release/brey/pyPoseidon)

This is a development project utilising multiple solvers (currently DELFT3D & SCHISM) for simulating sea level height (currently only storm surge). The purpose is to create a simple, portable and transparent way of setting up, running and analysing hydrodynamic computations through python scripts and Jupyter Notebooks (http://jupyter.org). See Notebooks in examples/ for relevant prototypes.

## Instalation


`conda install -c gbrey pyPoseidon`


### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. See Wiki for more details.

SCHISM needs to be compiled for your system. You can download it from  http://columbia.vims.edu/schism/tags/. See http://ccrm.vims.edu/schismweb/ for more info.


You can also install the solvers easily with conda  

`conda install -c gbrey pschism delft3d4`


## Tests

There are several sets of tests. You can run pyPoseidon unitests with 

`pytest --pyargs pyPoseidon`

In order to test also the solver integration use 

`python -m pyPoseidon.tests --runschism`

or

`python -m pyPoseidon.tests --rundelft`

if you are using a local installation of the solvers please specify the PATH to the executables in your system such as 

`export D3D = '/path_to_folder_bin/lnx64/flow2d3d/'`
`export LD3D = '/path_to_folder_bin/lnx64/flow2d3d/'`
`export SCHISM = '/path_to_schism_executable'`


## License
* The project is released under the EUPL v1.2 license. 
