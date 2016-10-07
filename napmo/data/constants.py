# file: constants.py
# nAPMO package
# Copyright (c) 2015,  Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

"""
UnitsDatabase
=============

This database contains factors to covert between Bohr and Angstroms.
And some other physical constants:

* ``BOHR_TO_ANGSTROM = 0.52917721092``
* ``ANGSTROM_TO_BOHR = 1.889725989``
* ``ELECTRON_CHARGE = -1.0``
* ``ELECTRON_MASS = 1.0``
* ``PROTON_CHARGE = 1.0``
* ``PROTON_MASS = 1836.15267247``
* ``NEUTRON_CHARGE = 0.0``
* ``NEUTRON_MASS = 1838.6836605``
* ``SPIN_UP_ELECTRON = 0.5``
* ``SPIN_DOWN_ELECTRON = -0.5``
* ``SPIN_ELECTRON = 0.5``
* ``SPIN_FERMION = 0.5``
* ``SPIN_BOSON = 1.0``

"""

BOHR_TO_ANGSTROM = 0.52917721092
ANGSTROM_TO_BOHR = 1.889726133921252

ELECTRON_CHARGE = -1.0
ELECTRON_MASS = 1.0
PROTON_CHARGE = 1.0
PROTON_MASS = 1836.15267247
NEUTRON_CHARGE = 0.0
NEUTRON_MASS = 1838.6836605

SPIN_UP_ELECTRON = 0.5
SPIN_DOWN_ELECTRON = -0.5
SPIN_ELECTRON = 0.5

# n * 0.5 n = 1, 2, 3, ..., n (U.A.)
SPIN_FERMION = 0.5

# n = 0, 1, 2, 3, ..., n (U.A.)
SPIN_BOSON = 1.0
