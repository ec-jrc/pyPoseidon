STORM SURGE ESTIMATION
==============================

This is a development project utilizing (for now) the DELFT3D open source code for considering storm surge. The purpose is to create, run and analyze a storm surge computation based on Python scripts and thus avoid the Matlab based GUI that Deltares uses. The procedure is outlined in Jupyter notebooks in order to have clarity and transparency.

## Instalation

The best way is with conda, see here (https://conda.io/docs/intro.html).

Then create a dedicated environment with

* conda env create -f environment.yml 

You should have a environment named Poseidon (check it with 'conda env list')

Activate the environment 

* source activate Poseidon

For using the package through jupyter notebooks

* conda config --add channels conda-forge
* conda install ipykernel
* python -m ipykernel install --name Poseidon

Then choose the Poseidon kernel on your Jupyter notebook. 


Finally, install Poseidon the usual way

* python setup.py install

### Prerequisities

DELFT3D needs to be compiled for your system. You can download it from http://oss.deltares.nl/web/delft3d/source-code. Follow the instruction therein.

For using the animation capabilities as outlined in the Notebooks, FFmpeg (https://ffmpeg.org/) library needs to be installed. 

## Tests

No tests are available at the moment.

## Authors

* **George Breyiannis** 


## Acknowledgments

Modified versions of modules from the project BOUT++ (https://github.com/boutproject/BOUT-dev) and OpenEarthTools (https://publicwiki.deltares.nl/display/OET/OpenEarth) are incorporated into Poseidon.

* All the people that teach me stuff.  

## License
* The project is released under the EUPL 1.1 license. 

