#!/usr/bin/env python3
# file: Slater_density.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division

import numpy as np
import os
import time

from copy import deepcopy

from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *

"""
Calcualtion of :math:`\int \\rho(\\bf r)` for several diatomic molecules.
"""


def init_system(element, distance, basis_kind, basis_file):
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom(element, [0.000000, 0.000000, distance/2.0], basis_kind=basis_kind, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance/2.0], basis_kind=basis_kind, basis_file=basis_file)

    # Functional definition (for Python)
    def rho(coord, molecule):
        basis = molecule.get_basis_set('e-')
        occupation = molecule.n_occupation('e-')
        bvalue = basis.compute(coord)
        output = 0.0

        for k in range(occupation):
            output += bvalue[k] * bvalue[k]

        return output * 2

    # Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')

    return molecule, exact, rho


if __name__ == '__main__':

    # Grid definition
    angularPoints = 110
    radialPoints = 50
    grid = BeckeGrid(radialPoints, angularPoints)
    grid.show()

    # Test for diatomic molecules for the following elements:
    elements = ['He', 'Li', 'O']
    distances = [2.75, 2.623, 1.206]
    basis_file = os.path.join(os.path.dirname(__file__), "STO.json")
    basis_kind = "STO"

    for (element, distance) in zip(elements, distances):
        molecule, exact, rho = init_system(element, distance, basis_kind, basis_file)

        # Calculate integral (Python Code)
        start_time = time.time()
        integral_p = grid.integrate(molecule, rho)
        elapsed_time_p = time.time() - start_time

        print("%4s %12.8f  %12.8f %12.7f" % (
            element+str(2), integral_p, np.abs(exact - integral_p), elapsed_time_p))

    grid.free()
