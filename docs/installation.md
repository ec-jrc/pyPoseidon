<style>body {text-align: justify}</style>

There are a number of ways to install `pyposeidon` and its dependencies.

### With `conda`

The easiest way to get started is to install the `pyposeidon` package from the `conda-forge`
channel:

```
conda install -c conda-forge pyposeidon
```

This will pull all necessary dependencies and is the recommended way to integrate `pyposeidon` to
your environment.


If you don't intend to use the visualisation modules of `pyposeidon` then, you may install the
`pyposeidon-base` package, instead, which is quite a bit lighter:

```
conda install -c conda-forge pyposeidon-base
```

### With `PyPI`

In order to install from `PyPI` please note that your system needs to have:

- `python>=3.8`
- `geos`
- `proj<8`
- `eccodes > 2.20`

Once you have those, you just need to run:

```
pip install pyposeidon
```

If you intend to use the visualisation modules, then you should install the necessary libraries with:

```
pip install pyposeidon --extras viz
```

!!! note

    Depending on which modules you intend to use you might need to install additional non-python
    dependencies (e.g. `jigsaw`, `gmesh` for mesh generation, `schism`, `delft3d` for numerical solving etc).


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

There are two flavours depending on the `mpi` option used. This can be explicitly selected with e.g.

```
conda install -c gbrey pschism=5.9=mpi_mpich_*
```

or

```
conda install -c gbrey pschism=5.9=mpi_openmpi_*
```

respectively.

Alternatively, you can download and compile them for your system:

- Delft3D can be downloaded from [deltares](http://oss.deltares.nl/web/delft3d/source-code).
- SCHISM can be downloaded from [github](https://github.com/schism-dev/schism).

## Mesh generation binaries

When installing `pyposeidon` via `PyPI` the mesh generation binaries are required and can be installed with  conda:

```
conda install -c conda-forge gmesh
conda install -c conda-forge jigsaw
```

Currently [Jigsaw](https://github.com/dengwirda/jigsaw) and [GMSH](http://gmsh.info) are supported.
