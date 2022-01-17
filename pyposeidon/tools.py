import logging
import os
import shlex
import pathlib
import subprocess
from collections.abc import Iterable
import time

import psutil
import xarray as xr
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

LAUNCH_SCHISM_TEMPLATE = """
#!/usr/bin/env bash
#
# launchschism.sh
#
# Launch schism using MPI

set -euo pipefail

mkdir -p outputs

exec mpirun {mpirun_flags} -N {ncores} {cmd}
""".strip()


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


def run(cmd: str, quiet: bool = False, check: bool = True) -> subprocess.CompletedProcess:
    t1 = time.perf_counter()
    proc = subprocess.run(shlex.split(cmd), check=False, capture_output=True, text=True)
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


def is_openmpi() -> bool:
    cmd = "mpirun --version"
    proc = run(cmd, quiet=True)
    return "Open MPI" in proc.stdout


def create_mpirun_script(
    target_dir: str,
    cmd: str,
    use_threads: bool = True,
    script_name: str = "launchSchism.sh",
    ncores: int = 0,
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
    if use_threads and is_openmpi():
        mpirun_flags = "--use-hwthread-cpus"
    else:
        mpirun_flags = ""
    content = LAUNCH_SCHISM_TEMPLATE.format(
        mpirun_flags=mpirun_flags,
        ncores=ncores,
        cmd=cmd,
    )
    # Write to disk and make executable
    script_path = target_dir + "/" + script_name
    with open(script_path, "w") as fd:
        fd.write(content)
    make_executable(script_path)
    return script_path


def open_dataset(source: os.PathLike, **kwargs) -> xr.Dataset:
    """
    A wrapper around `xr.open_dataset` that supports multiple engines

    This function is needed because `xr.open_dataset()` is not yet supporting e.g. rasterio
    """
    logger.info("extracting dem from %s\n", source)
    if isinstance(source, pathlib.Path):
        source = source.as_posix()
    if source.lower().endswith("tif"):  # GeoTiff
        data_array = xr.open_rasterio(source, parse_coordinates=True, **kwargs)
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
    return isinstance(obj, Iterable)


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
