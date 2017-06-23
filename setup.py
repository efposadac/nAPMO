#!/usr/bin/env python3
# file: setup.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co
"""
Script for installation of the python module.

Info:
    The C extension is builded outside this script, be sure to make it before this one!
"""
from setuptools import find_packages, Extension, setup, Command
import os


def read(fname):
    """Include README file in the package."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if __name__ == '__main__':

    setup(
        name="napmo",
        version="0.1",
        description="Numerical Any Particle Molecular Orbital",
        author="Fernando Posada",
        author_email="eposada@sissa.it",
        url="http://efposadac.github.io/nAPMO/",
        packages=find_packages(),
        package_data={'': ['*.json', '*.dat', '*.so', 'data/basis/*']},
        long_description=read('README.rst')
    )

    exit(0)
