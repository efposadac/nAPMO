#!/usr/bin/env python
# file: Coulomb.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import time
import os

from napmo.utilities.lebedev import *
from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.grids.poisson_solver import poisson_solver


def rho_ab(a, b, coord):
    return a.compute(coord) * b.compute(coord)

if __name__ == '__main__':
    #
    # Molecule definition
    basis_file = os.path.join(os.path.dirname(__file__), "TEST.json")
    molecule = MolecularSystem()
    molecule.add_atom("H", [0.0, 0.0, 0.371], basis_kind="GTO", basis_file=basis_file)
    # molecule.add_atom("H", [0.0, 0.0, -0.371], basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # basis set
    basis = molecule.get('atoms')[-1].get('basis')
    a = basis.get('function')[0]
    b = basis.get('function')[0]

    # Build grid
    angularPoints = 110
    radialPoints = 100

    grid = BeckeGrid(molecule, radialPoints, angularPoints)

    # Calculate potential
    lmax = int(lebedev_get_order(angularPoints) / 2)
    rho = rho_ab(a, b, grid.points)

    U = poisson_solver(grid, rho, lmax)

    v = np.zeros(grid.size, dtype=np.float64)

    for i in range(grid.ncenter):
        grid.evaluate_decomposition(
            U[i][:], grid.origin[i, :], v)

    integral = grid.integrate(rho * v)

    assert np.allclose(integral, 0.77461)
