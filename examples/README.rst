**Numerical Any particle Molecular Orbital (nAPMO) package.**

:Author: Fernando Posada Correa, MHPC, 2015

This folder contains some examples using the code. The first one is a ipython notebook called ``nAPMO.ipynb`` which contains information step by step ob how the code works.

The second example can be found in the file ``Density.py``. This file contains examples for the calculation of :math:`\int \rho(\\bf r)` for diatomic molecules for elements from Z = 1, up to Z = 8, with exception of the Helium. The outcome of this test must be:
::
	NF:    14 H2  Int:    2.000 Error:   0.0001 Time:    4.2388203
	NF:    56 Li2 Int:    6.000 Error:   0.0003 Time:    9.9827564
	NF:    56 Be2 Int:    8.000 Error:   0.0002 Time:    9.9273388
	NF:    56 B2  Int:   10.000 Error:   0.0005 Time:   10.0260947
	NF:    56 C2  Int:   12.001 Error:   0.0013 Time:   10.1134686
	NF:    56 N2  Int:   14.001 Error:   0.0009 Time:    9.9720314
	NF:    56 O2  Int:   16.001 Error:   0.0007 Time:   10.0695314

NF is the number of functions on the basis-set employed. The grid used for this calculation was 100 radial x 110 angular. It can be seen that the performance is too slow, taking in to account that for each molecule the code is calculating only one integral!. However as a prototype it can be said that the algorithm works and gives the expected result.

All ``*.dens`` files are density matrices to perform the integration.