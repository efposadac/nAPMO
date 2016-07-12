#!/usr/bin/env python
# file: test_decomposition.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import os

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.lebedev import lebedev_get_order
from napmo.grids.becke import BeckeGrid
from napmo.grids.cubic_spline import CubicSpline


def test_decomposition_grid():
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0], basis_name='STO-3G')
    print(molecule)

    # Grid definition
    angularPoints = 110
    radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    for radialPoints in radialList:
        lmax = int(lebedev_get_order(angularPoints) / 2)
        lsize = (lmax + 1) ** 2

        # Grids
        grid = BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)
        atgrid = grid.atgrids[-1]

        # Basis and functional
        basis = molecule.get_basis('He')
        a = basis.get('function')[0]
        b = basis.get('function')[0]

        def p_ab(coord):
            return a.compute(coord) * b.compute(coord)

        rho = np.empty(grid.size)
        rho = p_ab(grid.points)

        p_lm = atgrid.spherical_expansion(lmax, rho)

        # Build rho again:
        rho_n = np.zeros(grid.size)
        p_lm_cs = []
        idx = 0

        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                aux = np.array(p_lm[:, idx])
                p_lm_cs.append(
                    CubicSpline(
                        aux,
                        rtransform=atgrid.radial_grid.rtransform))
                idx += 1

        grid.evaluate_decomposition(0, p_lm_cs, rho_n)
        assert np.allclose(rho_n, rho)


def test_decomposition_atomic_grid():
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0], basis_name='STO-3G')
    print(molecule)

    # Grid definition
    angularPoints = 110
    radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    for radialPoints in radialList:
        lmax = int(lebedev_get_order(angularPoints) / 2)
        lsize = (lmax + 1) ** 2

        # Grids
        grid = BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)
        atgrid = grid.atgrids[-1]

        # Basis and functional
        basis = molecule.get_basis('He')
        a = basis.get('function')[0]
        b = basis.get('function')[0]

        def p_ab(coord):
            return a.compute(coord) * b.compute(coord)

        rho = np.empty(grid.size)
        rho = p_ab(grid.points)

        p_lm = atgrid.spherical_expansion(lmax, rho)

        # Build rho again:
        rho_n = grid.atgrids[-1].evaluate_expansion(lmax, p_lm)
        assert np.allclose(rho_n, rho)

# test_decomposition_grid()
# test_decomposition_atomic_grid()