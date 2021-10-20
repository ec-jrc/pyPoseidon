import logging
import os
import shlex
import subprocess
import time

import psutil

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
