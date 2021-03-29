pyPoseidon
==========

<style>body {text-align: justify}</style>

Hydrodynamic simulations within the context of geoflows are very cumbersome. The process can be split into 3 stages. 

The pre-processing stage where the necessary input is prepared. This includes input data conversion, parameter setting and grid generation. The latter usually requires a manual process using dedicated software.

The execution stage which invokes a numerical solver, usually written in fortran or C, that handles the numerics. 

Finally the post-processing stage where the output of the simulation is analyzed and/or visualized in order to gain insight and/or providing products & services.

The main purpose of pyPoseidon is to integrate the above stages. This is done by using python to handle the first and third step while wrapping around external code such as the numerical solvers, essentially treating them as plugins. This way the whole workflow is abstracted to the user. When python API is available, the integration is more straightforward. 

pyPoseidon strives to provide:

- *Reproducibility* through extensive logging, simulation signature, etc.

- *Transparency* by being an open source, well documented python package.

- *Expandability* by virtue of it's design and the ability to incorporate additional solvers.

- *Portability* as a feature with conda & python.

- *Scalability* via providing a consistent framework from a portable computer to HPC.

- *Interoperability* by depending on an ever expanding set of popular python structures like pandas, xarray etc.


Currently 2 numerical solvers are supported. 

- *Delft3D* which uses a structured grid. For more info see http://oss.deltares.nl/web/delft3d/source-code.

- *SCHISM* which relies on an unstructured grid. For more info see http://ccrm.vims.edu/schismweb/. 

# Documentation

## [Installation](10_installation.md) 

## [Structure](20_structure.md)

## [User Guide](20_usage.md)

## [API](50_API.md)

## Origin

* This project was initiated at the European Commission's Joint Research Centre within the E.1 unit during research into storm surge modelling.

## License
* The project is released under the EUPL 1.2 license.

