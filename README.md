#Numerical Any particle Molecular Orbital (nAPMO) package.#

* nAPMO is a numerical implementation of APMO - NEO - OMNE method.

The basic idea is to implement some of the numerical analysis tools required to perform a  numerical Hartree-Fock calculation as first step to implement the numerical OMNE approach.

* Version 0.0

### Compilation ###

At this moment there is a module in FORTRAN code. I use f2py to interface it with Python. To perform compilation just run:

* ``make``
* ``make clean`` is also supported.

One example using the code is provided in iphyton notebook format. (nAPMO-test.ipynb)

Any suggestion and help is more than welcome and appreciated. 

eposada@sissa.it