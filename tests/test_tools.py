from __future__ import annotations

import logging
import os
import pathlib
import shutil

import bottleneck as bn
import numcodecs.abc
import numpy as np
import pandas as pd
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
        ("schism develop", "develop"),
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
    cmd = "schism"
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
        assert f"/bin/{cmd} {scribes}" in cmd_line


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


def test_calc_quantization_params() -> None:
    params = tools.calc_quantization_params(0, 10, dtype=np.int8)
    assert pytest.approx(params["scale_factor"], abs=1e-3) == 0.039
    assert pytest.approx(params["add_offset"], abs=1e-3) == 5.039
    assert params["missing_value"] == 127
    assert params["dtype"] == np.int8
    #
    params = tools.calc_quantization_params(0, 10, dtype=np.int16)
    assert pytest.approx(params["scale_factor"], abs=1e-4) == 0.0001
    assert pytest.approx(params["add_offset"], abs=1e-4) == 5.0001
    assert params["missing_value"] == 32767
    assert params["dtype"] == np.int16
    #
    params = tools.calc_quantization_params(0, 10, dtype=np.int32)
    assert pytest.approx(params["scale_factor"], abs=1e-8) == 0
    assert pytest.approx(params["add_offset"], abs=1e-8) == 5
    assert params["missing_value"] == 2147483647
    assert params["dtype"] == np.int32


def test_quantization_roundtrip_with_nans() -> None:
    original = np.array([0, 4.9, np.nan, 10, np.nan])
    params = tools.calc_quantization_params(bn.nanmin(original), bn.nanmax(original), dtype=np.int8)
    expected = np.array([-128, -4, 127, 126, 127], dtype=np.int8)
    quantized = tools.quantize(original, **params)
    assert np.allclose(expected, quantized)
    dequantized = tools.dequantize(quantized, **params)
    assert np.allclose(original, dequantized, atol=1e-1, equal_nan=True)


def test_zarr_encoding_compressor_default():
    ds = pd.DataFrame({"a": [1, 5, 10, np.nan]}).to_xarray()
    default_compressor = numcodecs.Zstd(level=3)
    encoding = tools.get_zarr_encoding(ds)
    expected = {
        "a": {"compressor": default_compressor},
        "index": {"compressor": default_compressor},
    }
    assert encoding == expected


@pytest.mark.parametrize("compressor", [None, pytest.param(numcodecs.Blosc(cname="zstd", clevel=1), id="blosc")])
def test_zarr_encoding_compressor_custom(compressor):
    ds = pd.DataFrame({"a": [1, 5, 10, np.nan]}).to_xarray()
    encoding = tools.get_zarr_encoding(ds, compressor=compressor)
    expected = {
        "a": {"compressor": compressor},
        "index": {"compressor": compressor},
    }
    assert encoding == expected


def test_zarr_encoding_quantized_vars():
    ds = pd.DataFrame({"a": [1, 5, 10, np.nan]}).to_xarray()
    default_compressor = numcodecs.Zstd(level=3)
    encoding = tools.get_zarr_encoding(ds, quantized_vars={"a": np.int8})
    expected = {
        "a": {
            "compressor": default_compressor,
            "scale_factor": 0.03543307086614173,
            "add_offset": 5.535433070866142,
            "dtype": np.int8,
            "missing_value": 127,
        },
        "index": {"compressor": default_compressor},
    }
    assert encoding == expected


def test_netcdfr_encoding_quantized_vars():
    ds = pd.DataFrame({"a": [1, 5, 10, np.nan]}).to_xarray()
    ds = ds.chunk(index=1)
    encoding = tools.get_netcdf_encoding(ds, quantized_vars={"a": np.int8})
    expected = {
        "a": {
            "zlib": True,
            "complevel": 1,
            "chunksizes": [1],
            "scale_factor": 0.03543307086614173,
            "add_offset": 5.535433070866142,
            "dtype": np.int8,
            "missing_value": 127,
        },
        "index": {"zlib": True, "complevel": 1},
    }
    assert encoding == expected
