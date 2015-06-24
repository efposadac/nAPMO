#!/usr/bin/env python
# file: setup.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
"""
Script for installation of the python module.

It also builds the c library
"""
import os
from setuptools import find_packages, Extension, setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

napmo_module = Extension(
                      'libnapmo',
                      sources=['src/becke_grid.c', 'src/lebedev.c', 'src/gauss_chebyshev.c'],
                      extra_compile_args=['-O3', '-ffast-math', '-D _OMP', '-fopenmp'],
                      extra_link_args=['-lgomp']
                      )

setup(name="napmo",
      version="0.1",
      description="Numerical Any Particle Molecular Orbital",
      author="Fernando Posada",
      author_email="eposada@sissa.it",
      url="http://efposadac.github.io/nAPMO/",
      packages=find_packages(),
      package_data={'': ['*.json', '*.dat']},
      long_description=read('README.rst'),
      ext_modules=[napmo_module])
