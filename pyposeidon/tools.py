# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.
from __future__ import annotations

import itertools
import logging
import os
import pathlib
import re
import shlex
import shutil
import subprocess
import sys
import time
import typing as T
from collections import abc
from typing import TypeVar

import bottleneck as bn
import cartopy.feature
import colorlog
import geopandas as gpd
import jinja2
import numcodecs
import numpy as np
import numpy.typing as npt
import psutil
import rioxarray
import xarray as xr

from pyposeidon.utils.get_value import get_value


_PLAIN_FORMATTER = {
    "fmt": "%(asctime)s %(levelname)-8s %(name)s %(funcName)s:%(lineno)s %(message)s",
}
_COLORED_FORMATTER = {
    "fmt": "%(thin)s%(asctime)s%(reset)s %(log_color)s%(levelname)-8s%(reset)s %(thin)s%(name)s %(funcName)s:%(lineno)s%(reset)s %(message_log_color)s%(message)s",
    "secondary_log_colors": {
        "message": {"INFO": "bold_green", "WARNING": "bold_yellow", "ERROR": "bold_red", "CRITICAL": "bold_red"}
    },
}


logger = logging.getLogger(__name__)


SCHISM_VERSION_PATTERN = re.compile(r"schism v(\d+\.\d+\.\d+)\w*")


LAUNCH_SCHISM_TEMPLATE = """
#!/usr/bin/env bash
#
# launchschism.sh
#
# Launch schism using MPI

set -euo pipefail

if [[ ! -x "$(command -v realpath)" ]]; then
    realpath() {
        python -c "import os; print(os.path.realpath('$1'))"
    }
fi

root_dir="$(dirname "$(realpath "$0")")"

cd "${root_dir}"
mkdir -p outputs

# Schism sometimes throws an error mentioning ABORT or MPI_ABORT while it returns a status code of 0.
# In order to circumvent this we need to cacture the output and explicitly check the contents for ABORT.
schism_output=$(
    {{ mpirun_path }} {{ mpirun_flags }} -N {{ ncores }} {{ cmd }} {{ scribes }} 2>&1
)

echo "${schism_output}"

if [[ "${schism_output}" == *"ABORT"* ]]; then
  echo 'schism seems to have failed. Please check `err.log` and `run.log`. Exiting.'
  exit 111;
fi
""".strip()


LAUNCH_D3D_TEMPLATE = """
#!/bin/bash

export argfile=$1

if [ $# -gt 1 ]; then
       D3D="$2"
       LD3D="$3"
fi

#if [ -z $D3D ];then
#   echo 'no executable'
#   exit 1
#fi

    #
    # Set the directory containing delftflow.exe here
    #
exedir=$D3D/bin
libdir=$LD3D/lib

    #
    # No adaptions needed below
    #

    # Set some (environment) parameters
export LD_LIBRARY_PATH=$exedir:$libdir:$LD_LIBRARY_PATH
export PATH=$exedir:$PATH

    # Run
"$(which mpiexec)" {{ mpirun_flags }} -np {{ ncores }} {{ cmd }} $argfile
""".strip()


def setup_logging(
    min_level: int = logging.DEBUG, color: bool = True, log_file: os.PathLike[str] | str | None = "pyposeidon.log"
) -> None:
    # The purpose is to have a function that will allow us to easily setup some pyposeidon logging
    # and that will allow us to also dynamically change the log levels
    # The root logger and the other loggers should not be affected.

    # formatters
    colored_formatter = colorlog.ColoredFormatter(**_COLORED_FORMATTER)
    plain_formatter = logging.Formatter(**_PLAIN_FORMATTER)

    # Get handle of the root pyposeidon logger
    # we need to removing the existing handlers because if we call the function more than once
    # then we will end up with duplicated handlers (i.e. duplicated logs)
    logger = logging.getLogger("pyposeidon")
    logger.setLevel(min_level)
    logger.propagate = False
    logger.handlers.clear()

    # Setup the console
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG)
    if color:
        console_handler.setFormatter(colored_formatter)
    else:
        console_handler.setFormatter(plain_formatter)
    logger.addHandler(console_handler)

    # Setup the file handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(plain_formatter)
        logger.addHandler(file_handler)


def parse_schism_version(version_output: str) -> str:
    if "schism develop" in version_output:
        version = "develop"
    else:
        try:
            version_line = version_output.strip().splitlines()[0]
            version = SCHISM_VERSION_PATTERN.match(version_line).group(1)
        except Exception as exc:
            raise ValueError(f"Failed to parse version from:\n {version_output}") from exc
    return version


def get_schism_version() -> str:
    cmd = "schism -v"
    proc = subprocess.run(
        shlex.split(cmd),
        check=True,
        capture_output=True,
        text=True,
    )
    version = parse_schism_version(proc.stdout)
    return version


def get_solver(solver_name: str):
    # Panos: We do the imports here, because there is some sort of cyclical imports issue
    # and no time to properly solve it
    if solver_name == "schism":
        from .schism import Schism

        solver = Schism
    elif solver_name == "d3d":
        from .d3d import d3d

        solver = d3d
    else:
        raise ValueError(f"Unknown solver_name: {solver_name}")
    return solver


# From: https://stackoverflow.com/a/30463972/592289
def make_executable(path):
    mode = os.stat(path).st_mode
    mode |= (mode & 0o444) >> 2  # copy R bits to X
    os.chmod(path, mode)


def run(cmd: str, quiet: bool = False, check: bool = True, **kwargs) -> subprocess.CompletedProcess:
    t1 = time.perf_counter()
    proc = subprocess.run(shlex.split(cmd), check=False, capture_output=True, text=True, **kwargs)
    if not quiet:
        logger.info("Executed <%s> in %s seconds", cmd, time.perf_counter() - t1)
    if check and proc.returncode:
        # The process returned an error. Print stdout/stderr and raise
        logger.error("StdOut: %s", proc.stdout)
        logger.error("StdErr: %s", proc.stderr)
        proc.check_returncode()
    else:
        if not quiet:
            logger.info("StdOut: %s", proc.stdout)
            logger.info("StdErr: %s", proc.stderr)
    return proc


def is_mpirun_installed() -> bool:
    return bool(shutil.which("mpirun"))


def is_openmpi() -> bool:
    cmd = "mpirun --version"
    proc = run(cmd, quiet=True)
    return "Open MPI" in proc.stdout


def create_mpirun_script(
    template: str,
    script_name: str,
    target_dir: str,
    cmd: str,
    use_threads: bool = False,
    ncores: int = 0,
    scribes: int = -1,
) -> str:
    """
    Create a script for launching schism.

    - if `use_threads is True`, and the MPI implementation is `OpenMPI`, then the CPU threads are
      being used and `--use-hwthreaded_cpus` is passed to `mpirun`.
    - if `use_threads is False` or the MPI implementation is `mpich`, then only physical CPU cores
      are being used.
    """
    if ncores < 1:
        ncores = psutil.cpu_count(logical=use_threads)
    if not is_mpirun_installed():
        logger.error(f"mpirun is missing. Please install it and try again \n")
        return
    if use_threads and is_openmpi():
        mpirun_flags = "--use-hwthread-cpus"
    else:
        mpirun_flags = ""
    if scribes < 0:
        scribes = ""
    else:
        scribes = str(scribes)
    env = jinja2.Environment()
    template = env.from_string(template)
    content = template.render(
        mpirun_path=shutil.which("mpirun"),
        mpirun_flags=mpirun_flags,
        ncores=ncores,
        cmd=shutil.which(cmd),
        scribes=scribes,
    )
    # Write to disk and make executable
    script_path = os.path.join(target_dir, script_name)
    with open(script_path, "w") as fd:
        fd.write(content)
    make_executable(script_path)
    return script_path


def create_schism_mpirun_script(
    target_dir: str,
    cmd: str,
    use_threads: bool = False,
    script_name: str = "launchSchism.sh",
    template: str = LAUNCH_SCHISM_TEMPLATE,
    ncores: int = 0,
    scribes: int = -1,
) -> str:
    script_path = create_mpirun_script(
        target_dir=target_dir,
        cmd=cmd,
        use_threads=use_threads,
        script_name=script_name,
        ncores=ncores,
        scribes=scribes,
        template=template,
    )
    return script_path


def create_d3d_mpirun_script(
    target_dir: str,
    cmd: str = "d_hydro",
    use_threads: bool = False,
    script_name: str = "run_flow2d3d.sh",
    template: str = LAUNCH_D3D_TEMPLATE,
    ncores: int = 0,
) -> str:
    script_path = create_mpirun_script(
        target_dir=target_dir,
        cmd=cmd,
        use_threads=use_threads,
        script_name=script_name,
        ncores=ncores,
        template=template,
    )
    return script_path


def execute_schism_mpirun_script(cwd: str) -> None:
    """
    Execute launchSchism.sh and save stdout/stderr to disk
    """
    cmd = "./launchSchism.sh"
    with open(f"{cwd}/{cmd}", "r") as fd:
        contents = fd.read()
    proc = subprocess.run(
        shlex.split(cmd),
        check=False,
        capture_output=True,
        text=True,
        cwd=cwd,
    )

    with open(os.path.join(cwd, "err.log"), "w") as fd:
        fd.write(proc.stderr)
    with open(os.path.join(cwd, "run.log"), "w") as fd:
        fd.write(proc.stdout)

    if proc.returncode != 0:
        logger.error("schism failed to execute. Check the logs\n")
        proc.check_returncode()

    return proc


def open_dataset(source: os.PathLike, **kwargs) -> xr.Dataset:
    """
    A wrapper around `xr.open_dataset` that supports multiple engines
    """
    logger.info("extracting dem from %s\n", source)
    if isinstance(source, pathlib.Path):
        source = source.as_posix()
    if source.lower().endswith("tif"):  # GeoTiff
        data_array = rioxarray.open_rasterio(source, parse_coordinates=True, **kwargs)
        dataset = data_array.to_dataset(name="elevation").squeeze().reset_coords(drop=True)
    else:
        if source.lower().startswith("http"):  # URL
            engine = "pydap"
        else:
            engine = "netcdf4"
        logger.debug("Engine: %s", engine)
        dataset = xr.open_dataset(source, engine=engine, **kwargs)
    return dataset


def is_iterable(obj):
    return isinstance(obj, abc.Iterable)


def cast_path_to_str(path: os.PathLike) -> str:
    if isinstance(path, pathlib.Path):
        path = path.as_posix()
    return path


def flat_list(outer):
    return [item for inner in outer for item in inner]


def orient(nodes, tria, x="lon", y="lat"):
    # compute area
    ax = nodes.loc[tria.a, x].values
    bx = nodes.loc[tria.b, x].values
    cx = nodes.loc[tria.c, x].values
    ay = nodes.loc[tria.a, y].values
    by = nodes.loc[tria.b, y].values
    cy = nodes.loc[tria.c, y].values

    val = (by - ay) * (cx - bx) - (bx - ax) * (cy - by)

    tria["val"] = val

    tria["b"], tria["c"] = np.where(tria["val"] > 0, (tria["c"], tria["b"]), (tria["b"], tria["c"]))

    tria = tria.drop("val", axis=1)

    return nodes, tria


_RESOLUTIONS_TO_SCALES = {
    "HIGH": "10m",
    "High": "10m",
    "high": "10m",
    "H": "10m",
    "h": "10m",
    "MEDIUM": "50m",
    "Medium": "50m",
    "medium": "50m",
    "M": "50m",
    "m": "50m",
    "LOW": "110m",
    "Low": "110m",
    "low": "110m",
    "L": "110m",
    "l": "110m",
}


def get_coastlines(resolution: str, category="physical", name="land") -> gpd.GeoDataFrame:
    """
    Return a GeoDataFrame with the land mass at the specified ``resolution``.

    Wrapper of ``cartopy``'s ``NaturalEarthFeature``.
    """
    scale = _RESOLUTIONS_TO_SCALES[resolution]
    land = cartopy.feature.NaturalEarthFeature(
        category=category,
        name=name,
        scale=scale,
    )
    gdf = gpd.GeoDataFrame(geometry=list(land.geometries()))
    return gdf


# https://docs.python.org/3/library/itertools.html#itertools-recipes
# https://github.com/more-itertools/more-itertools/blob/2ff5943d76afa4591b5b4ae8cb4524578d365f67/more_itertools/recipes.pyi#L42-L48

_T = TypeVar("_T")
_U = TypeVar("_U")


def grouper(
    iterable: abc.Iterable[_T],
    n: int,
    *,
    incomplete: str = "fill",
    fillvalue: _U | None = None,
) -> abc.Iterator[tuple[_T | _U, ...]]:
    """Collect data into non-overlapping fixed-length chunks or blocks"""
    # grouper('ABCDEFG', 3, fillvalue='x') --> ABC DEF Gxx
    # grouper('ABCDEFG', 3, incomplete='strict') --> ABC DEF ValueError
    # grouper('ABCDEFG', 3, incomplete='ignore') --> ABC DEF
    args = [iter(iterable)] * n
    if incomplete == "fill":
        return itertools.zip_longest(*args, fillvalue=fillvalue)
    if incomplete == "strict":
        return zip(*args, strict=True)  # type: ignore[call-overload]
    if incomplete == "ignore":
        return zip(*args)
    else:
        raise ValueError("Expected fill, strict, or ignore")


def merge_netcdfs(paths: list[os.PathLike[str]], max_size: int = 5) -> xr.Dataset:
    # in order to keep memory consumption low, let's group the datasets
    # and merge them in batches
    datasets = [xr.open_dataset(path) for path in paths]
    while len(datasets) > max_size:
        datasets = [xr.merge(g for g in group if g) for group in grouper(datasets, max_size)]
    # Do the final merging
    ds = xr.merge(datasets)
    return ds


def resolve_schism_path(instance, kwargs) -> str:
    bin_path = os.environ.get("SCHISM", get_value(instance, kwargs, "epath", None))
    if bin_path is None:
        # ------------------------------------------------------------------------------
        logger.warning("Schism executable path (epath) not given -> using default \n")
        # ------------------------------------------------------------------------------
        bin_path = "schism"
    return bin_path


class QuantizationParams(T.TypedDict):
    scale_factor: float
    add_offset: float
    missing_value: int
    dtype: npt.DTypeLike


def calc_quantization_params(
    data_min: float,
    data_max: float,
    dtype: npt.DTypeLike,
) -> QuantizationParams:
    bits = np.iinfo(dtype).bits
    missing_value = (2**bits - 2) // 2
    scale_factor = (data_max - data_min) / (2**bits - 2)
    add_offset = data_min + 2 ** (bits - 1) * scale_factor
    return {
        "scale_factor": scale_factor,
        "add_offset": add_offset,
        "dtype": dtype,
        "missing_value": missing_value,
    }


def quantize(
    array: npt.NDArray[np.float_],
    *,
    add_offset: float,
    scale_factor: float,
    dtype: npt.DTypeLike,
    missing_value: int,
) -> npt.NDArray[np.int_]:
    nans = np.isnan(array)
    quantized: npt.NDArray[np.int_] = np.round((array - add_offset) / scale_factor, 0)
    quantized[nans] = missing_value
    return quantized.astype(dtype)


def dequantize(
    array: npt.NDArray[np.int_],
    *,
    add_offset: float,
    scale_factor: float,
    missing_value: int,
    dtype: npt.DTypeLike,
) -> npt.NDArray[np.float_]:
    array = bn.replace(array.astype(np.float64), missing_value, np.nan)
    dequantized = (array * scale_factor) + add_offset
    return dequantized


def update_or_add(outer: dict[str, T.Any], key: str, inner: abc.Mapping[str, T.Any]) -> None:
    if key in outer:
        outer[key].update(inner)
    else:
        outer[key] = inner


def get_zarr_encoding(
    ds: xr.Dataset,
    *,
    compressor: numcodecs.Blosc | None = numcodecs.Zstd(level=3),
    quantized_vars: abc.Mapping[str, npt.DTypeLike] | None = None,
) -> dict[str, dict[str, T.Any]]:
    encoding = ds.encoding.copy()
    for var in [str(var) for var in ds.variables]:
        update_or_add(encoding, var, {"compressor": compressor})
    if quantized_vars:
        for var, dtype in quantized_vars.items():
            params = calc_quantization_params(
                data_min=float(ds[var].min()),
                data_max=float(ds[var].max()),
                dtype=dtype,
            )
            update_or_add(encoding, var, params)
    return encoding


def get_netcdf_encoding(
    ds: xr.Dataset,
    *,
    quantized_vars: abc.Mapping[str, npt.DTypeLike] | None = None,
    compression_level: int = 1,
) -> dict[str, dict[str, T.Any]]:
    encoding = ds.encoding.copy()
    for var in [str(var) for var in ds.variables]:
        update_or_add(encoding, var, {"zlib": compression_level > 0, "complevel": compression_level})
        # Use chunks if they are defined
        if ds[var].chunks:
            chunksizes = [values[0] for values in ds[var].chunksizes.values()]
            update_or_add(encoding, var, {"chunksizes": chunksizes})
    if quantized_vars:
        for var, dtype in quantized_vars.items():
            params = calc_quantization_params(
                data_min=float(ds[var].min()),
                data_max=float(ds[var].max()),
                dtype=dtype,
            )
            update_or_add(encoding, var, params)
    return encoding
