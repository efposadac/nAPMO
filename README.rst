Numerical Any particle Molecular Orbital (NAPMO)
================================================

:Author: Fernando Posada Correa, MHPC, UNAL, Temple, 2015

.. image:: https://travis-ci.com/efposadac/nAPMO.svg?token=HUrCr32Dap17ppyzhhdd&branch=hf_dev
    :target: https://travis-ci.com/efposadac/nAPMO


NAPMO is a numerical implementation of APMO - NEO - OMNE method.

This program implements the  basis-set free  Hartree-Fock OMNE approach.

* Edwin Fernando Posada
* eposada@unal.edu.co
* Version 0.8

**Prerequisites:**

* Robust C++ compiler with C++11 support
* Recent Libint_ library
* GSL library
* Scipy
* Python 2.7 or 3.x

**Compilation:**


To compile and install the code use (it will install the package in ``$USER/.local``):

::

	make $(FLAVOR)

``$(FLAVOR)`` could be one of ``SERIAL``, ``OMP`` or ``CUDA``. Default is ``OMP``

To check the code run:

::

	nosetests --with-coverage --with-doctest

**Notes:**

This version contains:

* Multicenter molecular integrator based on Becke's paper
* Poisson solver based on Becke's strategy
* Analytical single or multi species Hartree-Fock solver
* Numerical single or multi species Hartree-Fock solver
* The C interface has been parallelized  with OMP and "some" CUDA. That code should be enabled in compilation time.

Any suggestion and help is more than welcome and appreciated. 

The code documentation can be found at http://efposadac.github.io/nAPMO/

fernando.posada@temple.edu


.. _libint: https://github.com/evaleev/libint