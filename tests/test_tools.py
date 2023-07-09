import logging
import os
import pathlib
import shutil

import psutil
import pytest
import xarray as xr

from pyposeidon import tools

from . import DATA_DIR


@pytest.mark.parametrize(
    "version_output,expected_version",
    [
        ("schism v5.9.0mod", "5.9.0"),
        ("schism v5.9.0", "5.9.0"),
        ("schism v5.10.1", "5.10.1"),
    ],
)
def test_parse_schism_version(version_output, expected_version):
    assert tools.parse_schism_version(version_output) == expected_version


@pytest.mark.skipif(not shutil.which("schism"), reason="requires schism binary")
def test_get_schism_version():
    version = tools.get_schism_version()
    assert isinstance(version, str)
    assert len(version) > 0


@pytest.mark.skipif(not shutil.which("mpirun"), reason="requires MPI backend")
@pytest.mark.parametrize("use_threads", [True, False])
@pytest.mark.parametrize("ncores", [0, 1, 2, 4, 44])
@pytest.mark.parametrize("scribes", [-1, 1])
def test_create_schism_mpirun_script(tmp_path, use_threads, ncores, scribes):
    target_dir = tmp_path.as_posix()
    cmd = "/bin/schism"
    script_name = "launchSchism.sh"
    script_path = tools.create_schism_mpirun_script(
        target_dir=target_dir,
        script_name=script_name,
        cmd=cmd,
        use_threads=use_threads,
        ncores=ncores,
        scribes=scribes,
    )
    assert script_path == (tmp_path / script_name).as_posix(), "script was not created"
    assert os.access(script_path, os.X_OK), "script is not executable"
    assert not (tmp_path / "outputs").exists(), "outputs subdirectory has not been created"
    content = pathlib.Path(script_path).read_text()
    for line in content.splitlines():
        if line.endswith("2>&1"):
            cmd_line = line
            break
    assert cmd_line
    assert cmd in cmd_line
    if ncores:
        assert f"-N {ncores}" in cmd_line
    else:
        assert f"-N {psutil.cpu_count(logical=use_threads)}" in cmd_line
    if scribes > 0:
        assert f"{cmd}) {scribes}" in cmd_line


@pytest.mark.skipif(not shutil.which("mpirun"), reason="requires MPI backend")
@pytest.mark.parametrize("ncores", [1, 2, 44])
def test_create_d3d_mpirun_script_ncores(tmp_path, ncores):
    target_dir = tmp_path.as_posix()
    cmd = "d_hydro"
    script_name = "run_flow2d3d.sh"
    script_path = tools.create_d3d_mpirun_script(
        target_dir=target_dir,
        script_name=script_name,
        cmd=cmd,
        ncores=ncores,
    )
    assert script_path == (tmp_path / script_name).as_posix(), "script was not created"
    assert os.access(script_path, os.X_OK), "script is not executable"
    assert not (tmp_path / "outputs").exists(), "outputs subdirectory has not been created"
    content = pathlib.Path(script_path).read_text()
    assert cmd in content, "the cmd is not being used"
    assert f"-np {ncores}" in content


@pytest.mark.parametrize(
    "source",
    [
        pytest.param(DATA_DIR / "dem.nc", id="netcdf - pathlib"),
        pytest.param(DATA_DIR / "dem.tif", id="geotiff - pathlib"),
        pytest.param((DATA_DIR / "dem.nc").as_posix(), id="netcdf- str"),
        pytest.param((DATA_DIR / "dem.tif").as_posix(), id="geotiff - str"),
        #        pytest.param("https://coastwatch.pfeg.noaa.gov/erddap/griddap/srtm30plus", id="URL"),
    ],
)
def test_open_dataset(source) -> None:
    ds = tools.open_dataset(source)
    assert isinstance(ds, xr.Dataset)


def test_setup_logging():
    logger = logging.getLogger("pyposeidon")
    # By default two handlers should be configured on pyposeidon
    assert len(logger.handlers) == 2
    # By default two log level should be DEBUG
    assert logger.level == logging.DEBUG
    # if we disable file handler, there should be just 1 handler
    tools.setup_logging(log_file=None)
    assert len(logger.handlers) == 1
    # we should be able to change the min_level
    tools.setup_logging(min_level=logging.INFO)
    assert logger.level == logging.INFO
