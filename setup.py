#!/usr/bin/env python3
# file: setup.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu
"""
Script for installation of the python module.

Info:
    The C extension is builded outside this script, be sure to make it before running this one!
"""
from setuptools import find_packages, setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.rst').read_text(encoding='utf-8')


setup(
    name="napmo",
    version="2.0.0",
    author="Fernando Posada",
    author_email="fernando.posada@temple.edu",
    description="Numerical Any Particle Molecular Orbital",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="http://efposadac.github.io/nAPMO/",
    license='GNU General Public License v3 (GPLv3)',
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    package_data={'': ['*.json', '*.dat', '*.so', 'data/basis/*']},
    python_requires=">=3.6",
    install_requires=['numpy', 'scipy', 'matplotlib'],
)
