# file: ntwobody.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


def compute_coulomb(grid, dens, lmax):
    """
    Calculates the Coulomb potential ``V`` for any density source by solving Poisson's equation
    following Becke procedure.

    This function takes into account the multi-center nature of the potential, so it decomposes
    the potential in spherical harmonics and build the potential using cubic spline interpolation
    see Becke's paper for details.

    Args:
        grid (BeckeGrid) : The molecular grid.
        dens (ndarray) : source of the density (grid.size).
        lmax (int) : Maximum ``l`` for spherical expansion.

    Return:
        J (ndarray) : The coulomb potential in the grid.

    See:
        napmo.poisson_solver
    """

    rho = dens.copy()
    with napmo.runtime.timeblock('Poisson solver'):
        U = napmo.poisson_solver(
            grid, rho, lmax)

    with napmo.runtime.timeblock('Interpolation'):
        J = np.zeros(grid.size, dtype=np.float64)
        for i in range(grid.ncenter):
            grid.evaluate_decomposition(i, U[i][:], J)

    return J
