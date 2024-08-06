# content of conftest.py

import pathlib

import pytest
import requests

from _pytest.mark import Mark

from . import DATA_DIR


EMPTY_MARK = Mark("", [], {})
RUNLAST_MARK = "runlast"


def should_run_last(item) -> bool:
    closest_marker = item.get_closest_marker(RUNLAST_MARK, default=EMPTY_MARK)
    return closest_marker is RUNLAST_MARK


def move_tests_to_the_end(items):
    start, end = [], []
    for item in items:
        if should_run_last(item):
            end.append(item)
        else:
            start.append(item)
    return start + end


def pytest_addoption(parser):
    parser.addoption("--runschism", action="store_true", default=False, help="run schism tests")
    parser.addoption("--rundelft", action="store_true", default=False, help="run delft tests")
    parser.addoption("--runslow", action="store_true", default=False, help="run slow tests")
    parser.addoption("--runviz", action="store_true", default=False, help="run viz tests")


def pytest_collection_modifyitems(config, items):
    should_run_schism = config.getoption("--runschism")
    should_run_delft = config.getoption("--rundelft")
    should_run_slow = config.getoption("--runslow")
    should_run_viz = config.getoption("--runviz")

    skip_schism = pytest.mark.skip(reason="need --runshism option to run")
    skip_delft = pytest.mark.skip(reason="need --rundelft option to run")
    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    skip_viz = pytest.mark.skip(reason="need --runviz option to run")

    run_last_marks = ("schism", "delft", "slow")
    for item in items:
        if any(mark in item.keywords for mark in run_last_marks):
            item.add_marker(RUNLAST_MARK)
        if "schism" in item.keywords and not should_run_schism:
            item.add_marker(skip_schism)
        if "delft" in item.keywords and not should_run_delft:
            item.add_marker(skip_delft)
        if "slow" in item.keywords and not should_run_slow:
            item.add_marker(skip_slow)
        if "viz" in item.keywords and not should_run_viz:
            item.add_marker(skip_viz)

    items[:] = move_tests_to_the_end(items)


def download_file_in_chunks(url: str, chunk_size: int = 1024):
    """Download file from `url` in chunks of size `chunk_size`."""
    response = requests.get(url, stream=True)
    response.raise_for_status()
    for chunk in response.iter_content(chunk_size=chunk_size):
        if chunk:  # filter out keep-alive new chunks
            yield chunk


def save_as(url: str, path: pathlib.Path) -> None:
    with path.open("wb") as fd:
        for chunk in download_file_in_chunks(url):
            fd.write(chunk)


TEST_DATA = {
    DATA_DIR
    / "ocean.parquet": "https://www.dropbox.com/scl/fi/msiru06qdh7cq99lsxicx/ocean.parquet?rlkey=ku7afiues47y91a15hdvtcybu&st=nayi0hc2&dl=1",
    DATA_DIR / "bl.zip": "https://www.dropbox.com/sh/9tfdl67sll1ax8c/AACntQeIavCzRfTZZ9Tp8uFda?dl=1",
    DATA_DIR / "dem.nc": "https://www.dropbox.com/s/l16crheqc9d89gy/dem.nc?dl=1",
    DATA_DIR / "dem.tif": "https://www.dropbox.com/s/dgdlnr2p5q66np7/dem.tif?dl=1",
    DATA_DIR / "era5.grib": "https://www.dropbox.com/s/mo8z8mv8k9kqj91/era5.grib?dl=1",
    DATA_DIR / "erai.grib": "https://www.dropbox.com/s/d56ctpo7nfptc1b/erai.grib?dl=1",
    DATA_DIR / "hgrid.gr3": "https://www.dropbox.com/s/73f1q13istbmacb/hgrid.gr3?dl=1",
    DATA_DIR / "meteo.nc": "https://www.dropbox.com/s/1kqt8732a284tbk/meteo.nc?dl=1",
    DATA_DIR / "uvp_2018100100.grib": "https://www.dropbox.com/s/jhf117jf9t2os8u/uvp_2018100100.grib?dl=1",
    DATA_DIR / "uvp_2018100112.grib": "https://www.dropbox.com/s/6y57hg119cpd6gx/uvp_2018100112.grib?dl=1",
    DATA_DIR / "uvp_2018100200.grib": "https://www.dropbox.com/s/7uejb2061k19zsf/uvp_2018100200.grib?dl=1",
    DATA_DIR / "uvp_2018100212.grib": "https://www.dropbox.com/s/nidbkobti83urzw/uvp_2018100212.grib?dl=1",
}


@pytest.hookimpl()
def pytest_sessionstart(session):
    for path, url in TEST_DATA.items():
        if not path.exists():
            print(f"Downloading: {url}")
            save_as(url, path)
