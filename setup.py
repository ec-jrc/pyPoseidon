#!/usr/bin/env python

from setuptools import setup, Extension
#from distutils.core import setup
#from distutils.extension import Extension
import numpy as np

def readme():
    with open('README.md') as f:
        return f.read()



setup(name='Poseidon',
      version='0.1',
      description='Storm Surge analysis tool',
      long_description=readme(),
      url='https://github.com/brey/Poseidon',
      author='George Breyiannis',
      author_email='gbreyiannis@gmail.com',
      license='EUPL',
      packages=['Poseidon.model','Poseidon.meteo','Poseidon.dem','Poseidon.tide','Poseidon.utils'],
      classifiers=[
          'Programming Language :: Python',
          'License :: OSI Approved :: EUPL',
          'Operating System :: OS Independent',
          'Development Status :: 4 - Beta',
          'Environment :: Other Environment',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering :: Atmospheric Science',
      ],
      ext_modules=[
          Extension("redtoreg", ["Poseidon/meteo/redtoreg.pyx"],
                    include_dirs=[np.get_include()]),
      ],
      package_data={'': ['misc/*']},
      zip_safe=False)

   
