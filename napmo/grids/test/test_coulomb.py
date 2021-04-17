#!/usr/bin/env python
# file: test_coulomb.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
import os

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.lebedev import *
from napmo.grids.becke import BeckeGrid
from napmo.grids.poisson_solver import poisson_solver


def rho_ab(a, b, coord):
    return a.compute(coord) * b.compute(coord)


def test_coulomb_H():
    #
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom("H", [0.0, 0.0, 0.0], basis_name='STO-3G')
    print(molecule)

    # basis set
    basis = molecule.get_basis('H')

    a = basis.get('cont')[0]
    b = basis.get('cont')[0]

    # Build grid
    angularPoints = 110
    radialPoints = 100

    grid = BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)

    # Calculate potential
    lmax = int(lebedev_get_order(angularPoints) / 2)
    rho = rho_ab(a, b, grid.points)

    U = poisson_solver(grid, rho, lmax)
    v = np.zeros(grid.size, dtype=np.float64)

    for i in range(grid.ncenter):
        grid.evaluate_decomposition(
            i, U[i][:], v)

    integral = grid.integrate(rho * v)
    # print(integral)
    assert np.allclose(integral, 0.77461)


if __name__ == '__main__':
    test_coulomb_H()
