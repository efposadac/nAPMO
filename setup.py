#!/usr/bin/env python3
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
from setuptools import find_packages, Extension, setup, Command
from distutils.command.build_ext import build_ext
from sys import argv, exit
import os
import subprocess


def read(fname):
    """Include README file in the package."""
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


def find_in_path(name, path):
    "Find a file in a search path"
    # adapted fom http://code.activestate.com/recipes/52224-find-a-file-given-a-search-path/
    for dir in path.split(os.pathsep):
        binpath = os.path.join(dir, name)
        if os.path.exists(binpath):
            return os.path.abspath(binpath)
    return None


def locate_cuda():
    """Locate the CUDA environment on the system

    Returns a dict with keys 'home', 'nvcc', 'include', and 'lib64'
    and values giving the absolute path to each directory.

    Starts by looking for the CUDAHOME env variable. If not found, everything
    is based on finding 'nvcc' in the PATH.
    """

    # first check if the CUDAHOME env variable is in use
    if 'CUDAHOME' in os.environ:
        home = os.environ['CUDAHOME']
        nvcc = os.path.join(home, 'bin', 'nvcc')
    else:
        # otherwise, search the PATH for NVCC
        nvcc = find_in_path('nvcc', os.environ['PATH'])
        if nvcc is None:
            raise EnvironmentError(
                'The nvcc binary could not be '
                'located in your $PATH. Either add it to your path, or set $CUDAHOME')

        home = os.path.dirname(os.path.dirname(nvcc))

    cudaconfig = {'home': home, 'nvcc': nvcc,
                  'include': os.path.join(home, 'include'),
                  'lib64': os.path.join(home, 'lib64')}

    for k, v in cudaconfig.items():
        if not os.path.exists(v):
            raise EnvironmentError('The CUDA %s path could not be located in %s' % (k, v))

    return cudaconfig


def customize_compiler_for_nvcc(self):
    """inject deep into distutils to customize how the dispatch
    to gcc/nvcc works.

    If you subclass UnixCCompiler, it's not trivial to get your subclass
    injected in, and still have the right customizations (i.e.
    distutils.sysconfig.customize_compiler) run on it. So instead of going
    the OO route, I have this. Note, it's kindof like a wierd functional
    subclassing going on."""

    # tell the compiler it can processes .cu
    self.src_extensions.append('.cu')

    # save references to the default compiler_so and _comple methods
    default_compiler_so = self.compiler_so
    super = self._compile

    # now redefine the _compile method. This gets executed for each
    # object but distutils doesn't have the ability to change compilers
    # based on source extension: we add it.
    def _compile(obj, src, ext, cc_args, extra_postargs, pp_opts):
        if os.path.splitext(src)[1] == '.cu':
            # use the cuda for .cu files
            self.set_executable('compiler_so', CUDA['nvcc'])
            # use only a subset of the extra_postargs, which are 1-1 translated
            # from the extra_compile_args in the Extension class
            postargs = extra_postargs['nvcc']
        else:
            postargs = extra_postargs['gcc']

        super(obj, src, ext, cc_args, postargs, pp_opts)
        # reset the default compiler_so, which we might have changed for cuda
        self.compiler_so = default_compiler_so

    # inject our redefined _compile method into the class
    self._compile = _compile


# run the customize_compiler
class custom_build_ext(build_ext):
    def build_extensions(self):
        customize_compiler_for_nvcc(self.compiler)
        build_ext.build_extensions(self)


class custom_build_type(Command):
    """Choose the type of build, OMP, CUDA or SERIAL"""

    description = "Defines what kind of extension you want to build, eg. CUDA, SERIAL or OMP"
    user_options = [('kind=', None, 'Kind of build; OMP, CUDA, SERIAL.')]

    def initialize_options(self):
        """init options"""
        self.kind = 'SERIAL'

    def finalize_options(self):
        pass

    def run(self):
        pass


if __name__ == '__main__':

    kind = 'SERIAL'
    if 'build_type' in argv:
        kind = argv[argv.index('build_type')+1].replace('=', ' ').split()[1]

    sources = ['src/becke_grid.c', 'src/lebedev.c', 'src/gauss_chebyshev.c']
    ccflags = ['-pg', '-g', '-Wall', '-O2', '-ffast-math', '-fPIC']
    ldflags = []
    libs_dir = []
    libs = []
    nvflags = []
    include = []

    if kind == 'OMP':
        ccflags += ['-D _OMP', '-fopenmp']
        ldflags += ['-lgomp']

    elif kind == 'CUDA':
        CUDA = locate_cuda()
        sources += ['src/gauss_chebyshev_cuda.cu']
        libs_dir += [CUDA['lib64']]
        libs += ['cudart', 'stdc++']
        include += [CUDA['include'], 'src']
        ccflags += ['-fopenmp', '-D _OMP', '-D _CUDA']
        ldflags += ['-lgomp']
        nvflags += ['-arch=sm_50']

        for flag in ccflags:
            nvflags.append('-Xcompiler')
            nvflags.append(flag)

    elif kind != 'SERIAL':
        print('build_type not allowed use CUDA, OMP or SERIAL')
        exit(1)

    napmo_module = Extension(
                          'libnapmo',
                          sources=sources,
                          library_dirs=libs_dir,
                          libraries=libs,
                          runtime_library_dirs=libs_dir,
                          extra_compile_args={'gcc': ccflags, 'nvcc': nvflags},
                          extra_link_args=ldflags,
                          include_dirs=include
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
          ext_modules=[napmo_module],
          cmdclass={'build_ext': custom_build_ext, 'build_type': custom_build_type})

    exit(0)
