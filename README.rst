Numerical Any particle Molecular Orbital (NAPMO)
================================================

:Author: Fernando Posada Correa, MHPC, UNAL, Temple U, 2015

.. image:: https://travis-ci.com/efposadac/nAPMO.svg?token=HUrCr32Dap17ppyzhhdd&branch=hf_dev
    :target: https://travis-ci.com/efposadac/nAPMO


NAPMO is a numerical implementation of APMO - NEO - OMNE method.

This program implements the  basis-set free  Hartree-Fock OMNE approach.

* Edwin Fernando Posada
* fernando.posada@temple.edu
* Version 1.1

**Prerequisites:**

* Robust C++ compiler with C++11 support
* Recent Libint_ library (v 2.6.0)
* GSL library
* Scipy
* Python 3.x
* pylibxc (v 4.3.4) provided by Libxc_

**Compilation:**

To compile and install the code use (it will install the package in ``$USER/.local``):

::

	make $(FLAVOR)

``$(FLAVOR)`` could be one of ``SERIAL``, ``OMP`` or ``CUDA``. Default is ``OMP``

To check the code run:

::

	nosetests --with-coverage --with-doctest --cover-package=napmo

**Notes:**

This version contains:

* GTO-based RHF/UHF and 
* GTO-based APMO-RHF/APMO-UHF
* Grid-based RHF/UHF 
* Grid-based APMO-RHF/APMO-UHF
* GTO-based RKS/UKS (LDA)
* GTO-based RKS-APMO (LDA/EPC17-2)

Any suggestion and help is more than welcome and appreciated. 

The code documentation can be found at http://efposadac.github.io/nAPMO/

fernando.posada@temple.edu

.. _libint: https://github.com/evaleev/libint
.. _libxc: https://www.tddft.org/programs/libxc/