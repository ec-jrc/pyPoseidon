import pathlib

from typing import List

import pyposeidon.meteo as pm
import pytest
import xarray as xr

from . import DATA_DIR


@pytest.fixture(scope="session")
def meteo_paths() -> List[pathlib.Path]:
    paths = list(sorted(DATA_DIR.glob("uvp_*")))
    return paths


@pytest.fixture(scope="session")
def meteo_datasets(meteo_paths) -> List[pm.meteo]:
    return [pm.meteo(meteo_source=path.as_posix(), meteo_engine="cfgrib").Dataset for path in meteo_paths]


def test_merge_strategy_last(meteo_paths, meteo_datasets):
    # In strategy "last" we want:
    # - the first 12 hours of all the datasets
    # - the rest of the hours of the last one
    expected = xr.concat(
        [
            *[ds.isel(time=slice(0, 12)) for ds in meteo_datasets],
            meteo_datasets[-1].isel(time=slice(12, None)),
        ],
        dim="time",
    )
    merged = pm.meteo(
        meteo_source=meteo_paths,
        meteo_engine="cfgrib",
        meteo_combine_by="nested",
        meteo_merge="last",
        meteo_xr_kwargs={"concat_dim": "step"},
    ).Dataset
    assert merged.equals(expected)


def test_merge_strategy_first(meteo_paths, meteo_datasets):
    # In strategy "first" we want:
    # - the first 13 hours of the first meteo
    # - hours 1-13 of all subsequent meteos
    # - hours 13-end of the last meteo
    expected = xr.concat(
        [
            meteo_datasets[0].isel(time=slice(0, 13)),
            *[ds.isel(time=slice(1, 13)) for ds in meteo_datasets[1:]],
            meteo_datasets[-1].isel(time=slice(13, None)),
        ],
        dim="time",
    )
    merged = pm.meteo(
        meteo_source=meteo_paths,
        meteo_engine="cfgrib",
        meteo_combine_by="nested",
        meteo_merge="first",
        meteo_xr_kwargs={"concat_dim": "step"},
    ).Dataset
    assert merged.equals(expected)
