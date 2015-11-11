#!/usr/bin/env python
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

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.utilities.density import *

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
    # molecule.show()

    # Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')

    # Get the density matrix (from a previous calculation)
    file_dens = os.path.join(os.path.dirname(__file__), element + '_dens.dat')

    return molecule, exact, file_dens

if __name__ == '__main__':

    # Grid definition
    angularPoints = 1202
    radialPoints = 1000

    # Test for diatomic molecules for the following elements:
    elements = ['H', 'O']
    distances = [0.742, 1.206]
    basis_name = "STO-3G"
    basis_file = os.path.join(os.path.dirname(__file__), "STO-3G.json")
    basis_kind = "GTO"

    # Header for results.
    print("System Int C         Error          Time C")

    for (element, distance) in zip(elements, distances):
        molecule, exact, file_dens = init_system(
            element, distance, basis_kind, basis_file)

        grid = BeckeGrid(molecule, radialPoints, angularPoints)
        # grid.show()

        basis = molecule.get_basis_as_cstruct('e-')

        # Calculate integral
        start_time = time.time()
        f = density_full_from_matrix_gto(file_dens, basis, grid.points)
        integral = grid.integrate(f)
        elapsed_time = time.time() - start_time

        # Print the results.
        print("%4s %12.8f  %12.8f  %12.7f" % (
            element + str(2), integral, np.abs(exact - integral), elapsed_time))
