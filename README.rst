**Numerical Any particle Molecular Orbital (nAPMO) package.**

:Author: Fernando Posada Correa, MHPC, 2015

nAPMO is a numerical implementation of APMO - NEO - OMNE method.

The basic idea is to implement the numerical tools required to perform a basis-set free Hartree-Fock calculation as a first step to implement the basis-set free OMNE approach.

* Edwin Fernando Posada
* eposada@sissa.it
* Version 0.1

Compilation
===========

To compile and install the code use (it will install the package in $USER/.local):

::

	make $(FLAVOR)

``$(FLAVOR)`` available are: SERIAL, OMP and CUDA.

To check the code run:

::

	nosetests --with-coverage --with-doctest

Notes
======

This version contains:

* Multicenter molecular integrator based on Becke's paper
* Poisson solver based on Becke's strategy
* Coulomb numerical integrals.
* Tools to manage a molecular system, basis-sets and others.
* All numerical methods use ctypes to accelerate the calculation.
* The C interface has been parallelized  with OMP and CUDA. That code can be enabled in compilation time.

Any suggestion and help is more than welcome and appreciated. 

Code documentation can be found at http://efposadac.github.io/nAPMO/

eposada@sissa.it

.. _examples: examples