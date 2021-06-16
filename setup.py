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
from setuptools import find_packages, setup, Extension
from setuptools.command.build_ext import build_ext
from distutils.version import LooseVersion

import subprocess
import pathlib
import os


class makeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = pathlib.Path(sourcedir).resolve()


class makeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['make', '--version'])
        except OSError:
            raise RuntimeError("Make must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        env = os.environ.copy()
        env['EXTDIR'] = extdir
        subprocess.check_call(['make', '-j4'], cwd=ext.sourcedir, env=env)


here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.rst').read_text(encoding='utf-8')


setup(
    name="napmo",
    version="2.0",
    author="Fernando Posada",
    author_email="fernando.posada@temple.edu",
    description="Numerical Any Particle Molecular Orbital",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="http://efposadac.github.io/nAPMO/",
    license='GNU General Public License v3 (GPLv3)',
    license_files="LICENSE.txt",
    ext_modules=[makeExtension('libnapmo', sourcedir='src')],
    cmdclass=dict(build_ext=makeBuild),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    scripts=['bin/napmo'],
    packages=find_packages(),
    package_data={'': ['*.json', '*.dat', '*.so', 'data/basis/*']},
    python_requires=">=3.6",
    install_requires=['numpy', 'scipy', 'matplotlib', 'pylibxc'],
)
