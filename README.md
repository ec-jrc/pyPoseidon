Framework for Hydrological simulations
=======================================

This is a development project utilising multiple solvers (currently DELFT3D, SCHISM) for considering storm surge. The purpose is to create a simple, portable and transparent way of setting up, running and analysing hydrodynamic computations through python scripts and Jupyter Notebooks (http://jupyter.org). See Notebooks in examples/ and Notebooks/ for relevant prototypes.

## Instalation

The best way is with conda, see here https://conda.io/docs/intro.html.

Download the Miniconda installer (python 3.6) for your system from https://conda.io/miniconda.html

Install Miniconda following the instructions here https://conda.io/docs/user-guide/install/index.html

Once your conda is installed (check with 'which conda' from your terminal), install Jupyter with 'conda install jupyter'.

Now you have a python3 base environment with jupyter support.

Then we create a dedicated environment for pyPoseidon. Navigate to where you have cloned the present repo and use

* conda env create -f environment.yml 

You should have an environment named pyPoseidon (check it with 'conda env list').

Activate the environment 

* source activate pyPoseidon

For using the environment in Jupyter notebooks use

* conda config --add channels conda-forge
* conda install ipykernel
* python -m ipykernel install --name pyPoseidon --user

Finally, install pyPoseidon the usual way

* python setup.py install

### Using Notebooks

From a new terminal window, start Jupyter with 

* jupyter notebook

Navigate to the pyPoseidon/examples & pyPoseidon/Notebook folders for using the Notebooks. 


### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. See Wiki for more details.
SCHISM needs to be compiled for your system. You can download it from  http://columbia.vims.edu/schism/tags/. See http://ccrm.vims.edu/schismweb/ for more info.

## Tests

No tests are available at the moment.


## License
* The project is released under the EUPL v1.2 license. 

