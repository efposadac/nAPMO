**Numerical Any particle Molecular Orbital (nAPMO) package.**

:Author: Fernando Posada Correa, MHPC, 2015

nAPMO is a numerical implementation of APMO - NEO - OMNE method.

The basic idea is to implement some of the numerical analysis tools required to perform a  numerical Hartree-Fock electronic calculation as first step to implement the numerical OMNE approach.

* Edwin Fernando Posada
* eposada@sissa.it
* Version 0.1

Compilation
===========

At this moment there is a module in FORTRAN code. I used f2py to interface it with Python. To perform compilation just run:
::
	make

To check the code go to the ``tests`` folder and run:
::
	nosetests --with-coverage --with-doctest

Notes
======

This version contains a complete and functional molecular multicenter integrator. At this moment it can be calculated overlap integrals over GTOs and STOs. As well the calculation of :math:`\int \rho(\bf r)`. Some examples are provided in the folder examples_.

Any suggestion and help is more than welcome and appreciated. 

eposada@sissa.it

.. _examples: examples