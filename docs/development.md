<style>body {text-align: justify}</style>

In order to develop `pyposeidon` you need:

1. Binary dependencies like e.g. `python`, `gdal`, `proj`. Additionally you will need solvers like
   `schism` and `d3d` and mesh generators like `jigsaw` and `gmsh`.
2. A bunch of runtime python dependencies like `numpy`, `xarray` and `matplotlib`.
3. Some python dependencies that are only useful during development, like `poetry` and `pytest`.

There are multiple ways to install all these. We are only going to suggest two of them:

## Install everything via conda

In the [locks](https://github.com/pmav99/pyPoseidon/tree/master/locks) directory there is a bunch of
[conda lock](https://github.com/conda-incubator/conda-lock) files. Choose one of the `dev` ones that
matches your Operating System and use it to create a new `conda` environment.

As an example, the following command will create a new conda environment named `pyposeidon_dev` and will install
*all* the `pyposeidon` dependencies.

```
env_name=pyposeidon_dev
conda create -n "${env_name}" --file conda-ubuntu-64-p38-openmpi-dev.lock
conda activate "${env_name}"
```

## Binary deps via conda + PyPI

Alternatively, you can choose one of the `binary` lock files and use it to create a base environment
with just the binary dependencies. After you create the environment you should use
[poetry](https://python-poetry.org/) to install the python dependencies.

```
env_name=pyposeidon_base
conda create -n "${env_name}" --file conda-ubuntu-64-p38-openmpi-binary.lock
conda activate "${env_name}"
poetry install
```

You are ready to go!


### Tests

You can use `pytest` to run the tests. Since the datasets required for testing `pyposeidon` are
substantial in size, they will be downloaded once you invoke `pytest` for the first time.

In order to execute the main test suite, it should be enough to use:

```
pytest
```

In order to test also the solver integration, use the following

```
pytest --runschism
```

 for testing `schism` or 

```
pytest --rundelft
```

for `delft3d` or both.

Some test are slow and can be invoked with 

```
pytest --runslow
```


If you are using a local installation for the solvers, please specify the `PATH` to the executables in your system such as:

```
export SCHISM = '/path_to_schism_executable'
```

or

```
export D3D = '/path_to_delft3d_bin_folder'
export LD3D = '/path_to_delft3d_lib_folder'
```
