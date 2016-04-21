#!/usr/bin/env python
# file: Slater_density.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division

import numpy as np
import os
import time

from copy import deepcopy

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid

"""
Calcualtion of :math:`\int \\rho(\\bf r)` for several diatomic molecules.
"""


def init_system(element, distance, basis_kind, basis_file):
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom(element, [0.000000, 0.000000, distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)

    # Functional definition (for Python)
    def rho(coord, molecule):
        basis = molecule.get_basis('e-')
        occupation = molecule.n_occupation('e-')
        bvalue = basis.compute(coord)

        output = 0.0
        for k in range(int(occupation)):
            output += bvalue[k] * bvalue[k]

        return output * 2

# Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')

    return molecule, exact, rho


if __name__ == '__main__':

    # Grid definition
    angularPoints = 194
    radialPoints = 100

    # Test for diatomic molecules for the following elements:
    elements = ['He', 'Li', 'O']
    distances = [2.75, 2.623, 1.206]
    basis_file = os.path.join(os.path.dirname(__file__), "STO.json")
    basis_kind = "STO"

    for (element, distance) in zip(elements, distances):
        molecule, exact, rho = init_system(
            element, distance, basis_kind, basis_file)

        grid = BeckeGrid(molecule, radialPoints, angularPoints)
        # grid.show()

        # Calculate integral
        start_time = time.time()
        f = rho(grid.points, molecule)
        integral = grid.integrate(f)
        elapsed_time = time.time() - start_time

        print("%4s %12.8f  %12.8f %12.7f" % (
            element + str(2), integral, np.abs(exact - integral), elapsed_time))
