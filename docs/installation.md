<style>body {text-align: justify}</style>

There are a number of ways to install `pyposeidon` and its dependencies

## `pyposeidon` installation

### `conda`

The easiest way to get started is to install the `pyposeidon` package from the `conda-forge`
channel:

```
conda install -c conda-forge pyposeidon
```

This will pull all necessary dependencies and is the recommended way to integrate `pyposeidon` to
your environment.

!!! note

    Due to an upstream issue, `pydap` needs to be installed manually. See
    `dependencies/python_deps10.yml` for info.

If you don't intend to use the visualization module of `pyposeidon` then, you may install the
`pyposeidon-base` package, instead, which is quite a bit lighter:

```
conda install -c conda-forge pyposeidon-base
```

### `PyPI`

In order to install from `PyPI` please note that your system needs to have:

- `python>=3.8`
- `geos`
- `proj<8`
- `eccodes > 2.20`

Once you have those, you just need to run:

```
pip install pyposeidon
```

If you intend to use the vizualization module, then you should install the necessary libraries with:

```
pip install pyposeidon --extras viz
```

!!! note

    Depending on which modules you intend to use you might need to install additional non-python
    dependencies (e.g. `jigsaw`, `gmesh` for mesh generation, `schism` for numerical solving etc).


## Solver integration

`pyposeidon` supports the following solvers:

- *[Delft3D](https://oss.deltares.nl/web/delft3d)*
- *[SCHISM](http://ccrm.vims.edu/schismweb/)*

!!! tip

    Support for `schism` is more mature.

The easiest way to install them is to use `conda`:

```
conda install -c gbrey pschism
conda install -c gbrey delft3d4
```

Alternatively, you can download them and compile them for your system:

- Delft3D can be downloaded from [deltares](http://oss.deltares.nl/web/delft3d/source-code).
- SCHISM can be downloaded from [github](https://github.com/schism-dev/schism).

## Mesh generation

`pyposeidon`'s `mesh` module can be used to generate localized or global meshes. To use the module you
must install one of the supported generators

```
conda install -c conda-forge gmesh
conda install -c conda-forge jigsaw
```
