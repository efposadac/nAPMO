#!/usr/bin/env python3
# file: Density.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division

import numpy as np
from copy import deepcopy

import os
import time

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
    # molecule.show()

    # Get the stack of atoms.
    atoms = molecule.get('atoms')

    # Build the C interface.
    system = 1  # CBinding(atoms)

    # Combine basis-set objects in one.
    basis = deepcopy(atoms[0].get('basis'))
    for i in range(1, len(atoms)):
        basis += atoms[i].get('basis')

    # Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')
    occupation = int(molecule.n_particles('e-') / 2)

    # Functional definition (for Python)
    def rho(coord, particle_stack, basis=basis):
        bvalue = basis.compute(coord)
        output = 0.0

        for k in range(occupation):
            output += bvalue[k] * bvalue[k]

        return output * 2

    return atoms, system, exact, rho


if __name__ == '__main__':

    # Grid definition
    angularPoints = 194
    radialPoints = 100
    grid = BeckeGrid(radialPoints, angularPoints)
    grid.show()

    # Test for diatomic molecules for the following elements:
    elements = ['He', 'Li', 'O']
    distances = [2.75, 2.623, 1.206]
    basis_name = "STO"
    basis_file = os.path.join(os.path.dirname(__file__), "STO.json")
    basis_kind = "STO"

    for (element, distance) in zip(elements, distances):
        atoms, system, exact, rho = init_system(element, distance, basis_kind, basis_file)

        # Calculate integral (Python Code)
        start_time = time.time()
        integral_p = grid.integrate(atoms, rho)
        elapsed_time_p = time.time() - start_time

        print("%4s %12.8f  %12.8f %12.7f" % (
            element+str(2), integral_p, np.abs(exact - integral_p), elapsed_time_p))

    grid.free()
