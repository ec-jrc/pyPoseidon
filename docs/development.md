<style>body {text-align: justify}</style>

In order to develop `pyposeidon` you need:

1. Binary dependencies like e.g. `python`, `gdal`, `proj`. Additionally you will need solvers like
   `schism` and `delft3d` and mesh generators like `jigsaw` and `gmesh`.
2. A bunch of runtime python dependencies like `numpy`, `xarray` and `matplotlib`.
3. Some python dependencies that are only useful during development, like `poetry` and `pytest`.

There are multiple ways to install all these. We are only going to suggest two of them:

## Install everything via conda

In the [locks](https://github.com/pmav99/pyPoseidon/tree/master/locks) directory there is a bunch of
[conda lock](https://github.com/conda-incubator/conda-lock) files. Choose one of the `dev` ones that
matches your Operating System and use it to create a new `conda` environment. E.g.

The following command will create a new conda environment named `pyposeidon_dev` and will install
*all* the `pyposeidon` dependencies **except** [pydap](https://github.com/pydap/pydap/), which must
be installed manually using `pip`.

!!! note

      The problem is that we need a recent `pydap` version and there has been no release for some
      time now. As soon as a new release is made, this should be resolved.

```
env_name=pyposeidon_dev
conda create -n "${env_name}" --file conda-ubuntu-64-p38-openmpi-dev.lock
conda activate "${env_name}"
pip install https://github.com/pydap/pydap/archive/3f6aa190c59e3bbc6c834377a61579b20275ff69.zip
```

## Binary deps via conda + PyPI

Alternatively, you can choose one of the `binary` lock files and use it to create a base environment
with just the binary dependencies. After you create the environment you should use
[poetry](https://python-poetry.org/) to install the python dependencies.

!!! tip

      In this case, you don't need to do anything special for
      [pydap](https://github.com/pydap/pydap/) since `poetry` handles installing packages straight
      from github.

```
env_name=pyposeidon_base
conda create -n "${env_name}" --file conda-ubuntu-64-p38-openmpi-binary.lock
conda activate "${env_name}"
poetry install
```

You are ready to go!


### Tests

You can use `pytest` to run the tests. Since the datasets required for testing `pyposeidon` are
substantial in size they will be downloaded once you invoke `pytest` for the first time.

In order to execute the main test suite, it should be enough to execute:

```
pytest
```

In order to test also the solver integration use the following:

```
pytest --runschism
```

If you are using a local installation of the solvers please specify the `PATH` to the executables in
your system such as:

```
export SCHISM = '/path_to_schism_executable'
```
