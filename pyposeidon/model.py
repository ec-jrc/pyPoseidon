"""
Main model module of pyposeidon. It controls the creation, execution & output of a complete simulation based on different hydrodynamic models

Currently supported : DELFT3D , SCHISM

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.

import os
import datetime
import numpy as np
import xml.dom.minidom as md
from shutil import copy2
import sys
import json
from collections import OrderedDict
import pandas as pd
import glob
from shutil import copyfile
import xarray as xr
import json

# local modules
from . import tools
import pyposeidon
import pyposeidon.mesh as pmesh
import pyposeidon.meteo as pmeteo
from pyposeidon.utils.get_value import get_value
import logging
import logging.config
import colorlog


def set(solver_name, atm=True, tide=False, **kwargs):
    """
    Construct a hydrodynamic model based on different solvers.
    !!! danger ""
        Due to a limitation of the Library rendering the docstrings, all arguments are marked
        as `required`, nevertheless they are all `Optional` except geometry.

    Args:

        solver_name str: Name of solver used, e.g. `d3d` or `schism`.
        atm bool: Flag for using meteo forcing.  Defaults to `True`.
        tide bool: Flag for using tidal configuration.  Defaults to `False`.

    """

    solver_class = tools.get_solver(solver_name=solver_name)
    instance = solver_class(atm=atm, tide=tide, **kwargs)
    return instance


def read(filename, **kwargs):
    end = filename.split(".")[-1]

    if end in ["txt", "json"]:
        with open(filename, "rb") as f:
            data = json.load(f)
            data = pd.json_normalize(data, max_level=0)
            info = data.to_dict(orient="records")[0]
    else:
        logger.error("Model file should be .txt or .json")
        sys.exit(0)

    return set(**info)
