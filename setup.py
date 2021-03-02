#!/usr/bin/env python

from setuptools import find_packages, setup
#from distutils.core import setup
#from distutils.extension import Extension
#import numpy as np

def readme():
    with open('README.md') as f:
        return f.read()

# copy structure from xarray/xarray-extras setup

DISTNAME = 'pyPoseidon'

ISRELEASED = False
VERSION = '0.4.1'
QUALIFIER = ''



LICENSE = 'EUPL 1.2'

AUTHOR = 'George Breyiannis'

AUTHOR_EMAIL = 'pyPoseidon@gmail.com'


CLASSIFIERS = [
    'License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)',
    'Operating System :: OS Independent',
    'Development Status :: 4 - Beta',
    'Environment :: Other Environment',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9'
        ]

DESCRIPTION = ''

LONG_DESCRIPTION = ''

PYTHON_REQUIRES = ''

INSTALL_REQUIRES = ''

SETUP_REQUIRES = ''


TESTS_REQUIRE = ['pytest >= 3.6']

URL = ''


setup(name=DISTNAME,
      version=VERSION,
      license=LICENSE,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      classifiers=CLASSIFIERS,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      python_requires=PYTHON_REQUIRES,
      install_requires=INSTALL_REQUIRES,
      setup_requires=SETUP_REQUIRES,
      tests_require=TESTS_REQUIRE,
      url=URL,
      packages=find_packages(),
      package_data={'pyPoseidon': ['misc/*','tests/data/*']})
