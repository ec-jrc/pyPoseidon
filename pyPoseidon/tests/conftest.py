# content of conftest.py

import pytest


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