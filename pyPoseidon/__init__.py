# -*- coding: utf-8 -*-

__author__ = 'George Breyiannis (breyiannis@gmail.com)'
__license__ = 'EUPL-1.2'
__version__ = '0.2.0'

__all__ = [
    'model',
    'meteo',
    'dem',
    'grid',
    'utils',
    'tide'
    'TEST_DIR',
]

import pathlib

from .model import *
from .d3d import d3d
from .schism import schism

# specify package paths
ROOT_DIR = pathlib.Path(__file__).resolve().parent
TEST_DIR = ROOT_DIR / "tests"
DATA_DIR = TEST_DIR / "data"
