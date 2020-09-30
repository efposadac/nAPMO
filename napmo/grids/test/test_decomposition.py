#!/usr/bin/env python
# file: test_decomposition.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

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
    molecule.add_atom("O", [0.0, 0.0, -1.0], basis_name='6-311G')
    print(molecule)

    basis = molecule.get_basis('e-')

    # Grid definition
    angularPoints = 110
    radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    for radialPoints in radialList:
        print("radial points", radialPoints)

        lmax = int(lebedev_get_order(angularPoints) / 2)
        lsize = (lmax + 1) ** 2

        # Grids
        grid = BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)

        rho = basis.compute(grid.points)
        rho = rho.sum(axis=1) * grid.becke_weights

        offset = 0
        p_lm_cs = []
        res = []
        for i in range(grid.ncenter):
            atgrid = grid.atgrids[i]
            p = rho[offset:offset + atgrid.size].copy()
            p_lm = atgrid.spherical_expansion(lmax, p)

            idx = 0
            for l in range(lmax + 1):
                for m in range(-l, l + 1):
                    aux = np.array(p_lm[:, idx])
                    p_lm_cs.append(CubicSpline(
                        aux, rtransform=atgrid.radial_grid.rtransform))
                    idx += 1
            offset += atgrid.size
            res.append(p_lm_cs)

        rho_n = np.zeros(grid.size)
        for i in range(grid.ncenter):
            grid.evaluate_decomposition(i, res[i][:], rho_n)
        assert np.allclose(rho_n, rho)

# test_decomposition_grid()
