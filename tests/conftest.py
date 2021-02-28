# content of conftest.py

import pathlib

import pytest
import requests

from . import DATA_DIR


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: mark test to run only when --slow is used"
    )
    config.addinivalue_line(
        "markers", "schism: mark test to run only when --runschism is used"
    )
    config.addinivalue_line(
        "markers", "delft: mark test to run only when --rundelft is used"
    )


def pytest_addoption(parser):
    parser.addoption(
        "--runschism", action="store_true", default=False, help="run schism tests"
    )
    parser.addoption(
        "--rundelft", action="store_true", default=False, help="run delft tests"
    )
    parser.addoption(
        "--slow", action="store_true", default=False, help="run slow tests"
    )


def pytest_collection_modifyitems(config, items):
    should_run_schism = config.getoption("--runschism")
    should_run_delft = config.getoption("--rundelft")
    should_run_slow = config.getoption("--slow")

    skip_schism = pytest.mark.skip(reason="need --runshism option to run")
    skip_delft = pytest.mark.skip(reason="need --rundelft option to run")
    skip_slow = pytest.mark.skip(reason="need --slow option to run")

    for item in items:
        if "schism" in item.keywords and not should_run_schism:
            item.add_marker(skip_schism)
        if "delft" in item.keywords and not should_run_delft:
            item.add_marker(skip_delft)
        if "slow" in item.keywords and not should_run_slow:
            item.add_marker(skip_slow)


def download_file_in_chunks(url: str, chunk_size: int = 1024):
    """ Download file from `url` in chunks of size `chunk_size`.  """
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
    DATA_DIR / "coast.cpg" : "https://www.dropbox.com/s/16841w9qpe07ulu/coast.cpg?dl=1",
    DATA_DIR / "coast.dbf" : "https://www.dropbox.com/s/m7ucjmg4l5b2gfl/coast.dbf?dl=1",
    DATA_DIR / "coast.prj" : "https://www.dropbox.com/s/dco2c8ifz5riahn/coast.prj?dl=1",
    DATA_DIR / "coast.shp" : "https://www.dropbox.com/s/dco2c8ifz5riahn/coast.prj?dl=1",
    DATA_DIR / "coast.shx" : "https://www.dropbox.com/s/695oekw7aebx73f/coast.shx?dl=1",
    DATA_DIR / "dem.nc" : "https://www.dropbox.com/s/l16crheqc9d89gy/dem.nc?dl=1",
    DATA_DIR / "era5.grib" : "https://www.dropbox.com/s/mo8z8mv8k9kqj91/era5.grib?dl=1",
    DATA_DIR / "erai.grib" : "https://www.dropbox.com/s/d56ctpo7nfptc1b/erai.grib?dl=1",
    DATA_DIR / "hgrid.gr3" : "https://www.dropbox.com/s/73f1q13istbmacb/hgrid.gr3?dl=1",
    DATA_DIR / "meteo.nc" : "https://www.dropbox.com/s/1kqt8732a284tbk/meteo.nc?dl=1",
    DATA_DIR / "param.in" : "https://www.dropbox.com/s/vf5mn7o7bv7yb5g/param.in?dl=1",
    DATA_DIR / "uvp_2018100100.grib" : "https://www.dropbox.com/s/jhf117jf9t2os8u/uvp_2018100100.grib?dl=1",
    DATA_DIR / "uvp_2018100112.grib" : "https://www.dropbox.com/s/6y57hg119cpd6gx/uvp_2018100112.grib?dl=1",
    DATA_DIR / "uvp_2018100200.grib" : "https://www.dropbox.com/s/7uejb2061k19zsf/uvp_2018100200.grib?dl=1",
    DATA_DIR / "uvp_2018100212.grib" : "https://www.dropbox.com/s/nidbkobti83urzw/uvp_2018100212.grib?dl=1",
    DATA_DIR / "rotated" / "E_JRC0000lf0000000020180101" : "https://www.dropbox.com/s/hs9wby5g939fqkn/E_JRC0000lf0000000020180101?dl=1",
    DATA_DIR / "rotated" / "E_JRC0000lf0001000020180101" : "https://www.dropbox.com/s/k3720x20hn9qw5b/E_JRC0000lf0001000020180101?dl=1",
    DATA_DIR / "rotated" / "E_JRC0000lf0002000020180101" : "https://www.dropbox.com/s/ta7yblif0n0vcsi/E_JRC0000lf0002000020180101?dl=1",
}


@pytest.hookimpl()
def pytest_sessionstart(session):
    for path, url in TEST_DATA.items():
        if not path.exists():
            print(f"Downloading: {url}")
            save_as(url, path)
