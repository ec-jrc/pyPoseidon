Framework for Hydrological simulations
=======================================

This is a development project utilizing (for now) the DELFT3D open source code for considering storm surge. The purpose is to create, run and analyze a storm surge computation based on Python scripts and thus avoid the Matlab based GUI that Deltares uses. The procedure is outlined in Jupyter notebooks in order to have clarity and transparency.

## Instalation

The best way is with conda, see here (https://conda.io/docs/intro.html).

Then create a dedicated environment with

* conda env create -f environment.yml 

You should have an environment named pyPoseidon (check it with 'conda env list').

Activate the environment 

* source activate pyPoseidon

For using the environment in jupyter notebooks

* conda config --add channels conda-forge
* conda install ipykernel
* python -m ipykernel install --name pyPoseidon

Then choose the pyPoseidon kernel on your Jupyter notebook. 

Finally, install pyPoseidon the usual way

* python setup.py install

### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. See Wiki for more details.

## Tests

No tests are available at the moment.

## Authors

* **George Breyiannis** 


## Acknowledgments

* All the people that teach me stuff.  

## License
* The project is released under the GPL v3 license. 

