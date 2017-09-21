#!/usr/bin/env python

from setuptools import setup, Extension
#from distutils.core import setup
#from distutils.extension import Extension
import numpy as np

def readme():
    with open('README.md') as f:
        return f.read()



setup(name='pyPoseidon',
      version='0.1',
      description='Storm Surge analysis tool',
      long_description=readme(),
      url='https://github.com/brey/pyPoseidon',
      author='George Breyiannis',
      author_email='gbreyiannis@gmail.com',
      license='GPLv3+',
      packages=['pyPoseidon', 'pyPoseidon.model','pyPoseidon.meteo','pyPoseidon.dem','pyPoseidon.utils','pyPoseidon.tide'],
      classifiers=[
          'Programming Language :: Python',
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Operating System :: OS Independent',
          'Development Status :: 4 - Beta',
          'Environment :: Other Environment',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      ext_modules=[
          Extension("redtoreg", ["pyPoseidon/meteo/redtoreg.pyx"],
                    include_dirs=[np.get_include()]),
      ],
      package_data={'': ['misc/*']})#,
#     zip_safe=False)

   
