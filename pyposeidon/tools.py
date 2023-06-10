# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import logging
import os
import shlex
import pathlib
import re
import subprocess
from collections.abc import Iterable
import time
import shutil

import psutil
import xarray as xr
import jinja2
import numpy as np
import rioxarray

logger = logging.getLogger(__name__)


SCHISM_VERSION_PATTERN = re.compile(r"schism v(\d+\.\d+\.\d+)\w*")


LAUNCH_SCHISM_TEMPLATE = """
#!/usr/bin/env bash
#
# launchschism.sh
#
# Launch schism using MPI

set -euo pipefail

root_dir="$(dirname "$(realpath "$0")")"

cd "${root_dir}"
mkdir -p outputs

exec $(which mpirun) {{ mpirun_flags }} -N {{ ncores }} $(which {{ cmd }}) {{ scribes }}
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


# TODO Handle master/develop version
def parse_schism_version(version_output: str) -> str:
    try:
        version_line = version_output.strip().splitlines()[0]
        version = SCHISM_VERSION_PATTERN.match(version_line).group(1)
        return version
    except Exception as exc:
        raise ValueError(f"Failed to parse version from:\n {version_output}") from exc


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
        mpirun_flags=mpirun_flags,
        ncores=ncores,
        cmd=cmd,
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
