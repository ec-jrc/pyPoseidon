# content of conftest.py

import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--with--solvers", action="store_true", default=False, help="run tests with solvers"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "solvers: mark test as need solvers to run")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--with--solvers"):
        # --with--solvers given in cli: do not skip solver tests
        return
    skip_solvers = pytest.mark.skip(reason="need --with--solvers option to run")
    for item in items:
        if "solvers" in item.keywords:
            item.add_marker(skip_solvers)
