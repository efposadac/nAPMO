# file: test_multicenter_integrator.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np
from copy import deepcopy
import os

from napmo.system.molecular_system import MolecularSystem
from napmo.system.cext import *
from napmo.solver.hf_solver import HF
from napmo.grids.becke import BeckeGrid


def test_multicenter_integrator():

    # Test for diatomic molecules of:
    elements = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O']
    distances = [0.742, 2.623, 2.427, 1.586, 1.268, 1.098, 1.206]
    basis_name = "6-31G"

    results = np.array([1.9999999595998637, 6.000002321070787,
                        7.999997958358362, 9.999998503574844,
                        11.999997997974441, 13.999997970177285,
                        16.00000891580461])
    count = 0

    for (element, distance) in zip(elements, distances):

        # print(element)

        # Molecule definition
        molecule = MolecularSystem()
        molecule.add_atom(
            element, [0.0, 0.0, distance / 2.0], basis_name=basis_name)
        molecule.add_atom(
            element, [0.0, 0.0, -distance / 2.0], basis_name=basis_name)

        # print(molecule)

        # Get the density matrix (from a previous calculation)
        hf = HF(molecule)
        energy = hf.compute(pprint=False)

        D = hf.PSI[-1].D

        # Grid definition
        angularPoints = 194
        radialPoints = 100
        grid = BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)
        # grid.show()

        basis = molecule.get_basis('e-')
        f = basis.compute(grid.points)
        f = np.array([aux.dot(D.dot(aux)) for aux in f])

        # Calculate integral
        integral = grid.integrate(f)
        print("INTEGRAL: ", integral, results[count])

        np.testing.assert_allclose(integral, results[count], rtol=10e-6)

        count += 1

# test_multicenter_integrator()
