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

    # Molecule definition
    basis_file = os.path.join(os.path.dirname(__file__), "TEST.json")
    molecule = MolecularSystem()
    molecule.add_atom("H", [0.0, 0.0, 0.0],
                      basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # basis set
    basis = molecule.get('atoms')[-1].get('basis')
    a = basis.get('function')[0]
    b = basis.get('function')[0]

    # Build grid
    angularList = [1202]
    radialList = [2000]

    for (angularPoints, radialPoints) in zip(angularList, radialList):
        grid = BeckeGrid(molecule, radialPoints, angularPoints)

        # Calculate potential
        start_time = time.time()
        lmax = int(lebedev_get_order(angularPoints) / 2)
        rho = rho_ab(a, b, grid.points)
        U = poisson_solver(grid, rho, lmax)

        # TODO: build molecular potential
        for i in range(grid.atgrids[-1].radial_grid.size):
            U[i, :] /= grid.atgrids[-1].radial_grid.points[i]

        v = grid.atgrids[-1].evaluate_expansion(lmax, U)

        integral = grid.integrate(rho * v)
        elapsed_time = time.time() - start_time

        print('$', radialPoints, '\\times', angularPoints,
              '$', '&', '$', integral, '$', '\\\\', elapsed_time)
