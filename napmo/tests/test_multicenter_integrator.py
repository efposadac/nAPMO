# file: test_multicenter_integrator.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np
from copy import deepcopy
import os

from napmo.system.molecular_system import MolecularSystem
from napmo.system.cext import *
from napmo.grids.becke import BeckeGrid
from napmo.utilities.density import *


def test_multicenter_integrator():

    # Test for diatomic molecules of:
    elements = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O']
    distances = [0.742, 2.623, 2.427, 1.586, 1.268, 1.098, 1.206]
    basis_name = "STO-3G"
    basis_file = os.path.join(
        os.path.dirname(__file__), "STO-3G.json")
    basis_kind = "GTO"

    results = np.array([1.9999999595998637, 6.000002321070787, 7.999997958358362,
                        9.999998503574844, 11.999997997974441, 13.999997970177285, 16.00000891580461])
    count = 0

    for (element, distance) in zip(elements, distances):

        print(element)

        # Molecule definition
        molecule = MolecularSystem()
        molecule.add_atom(element, [0.0, 0.0, distance / 2.0],
                          basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)
        molecule.add_atom(element, [0.0, 0.0, -distance / 2.0],
                          basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)

        # Get the density matrix (from a previous calculation)
        file_dens = os.path.join(os.path.dirname(
            __file__), element + '_dens.dat')

        # Grid definition
        angularPoints = 194
        radialPoints = 100
        grid = BeckeGrid(molecule, radialPoints, angularPoints)
        grid.show()

        basis = molecule.get_basis_as_cstruct('e-')

        # Calculate integral
        f = density_full_from_matrix_gto(file_dens, basis, grid.points)
        integral = grid.integrate(f)

        np.testing.assert_allclose(integral, results[count], rtol=10e-6)

        print(integral)

        count += 1

# test_multicenter_integrator()
