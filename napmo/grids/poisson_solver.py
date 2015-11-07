# file: poisson_solver.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.grids.radial import RadialGrid
from napmo.system.c_binding import napmo_library


def poisson_solver(grid, rho, lmax):
    offset = 0
    for i in range(grid.ncenter):
        atgrid = grid.atgrids[i]

        p = rho[offset:offset + atgrid.size]
        bw = grid.becke_weights[offset:offset + atgrid.size]

        # Spherical expansion
        sph_expansion = atgrid.spherical_expansion(lmax, p)

        # Boundary condition
        q = atgrid.integrate(p, bw)
        u_00 = 0.0
        if np.abs(q) > 1.0e-16:
            u_00 = np.sqrt(4.0 * np.pi * q)

        # Solve differential equation Ax = b
        U = np.zeros(sph_expansion.shape)
        b = np.zeros(atgrid.radial_grid.size + 2)
        A = np.zeros([atgrid.radial_grid.size + 2,
                      atgrid.radial_grid.size + 2])

        idx = 0
        pi4 = -4 * np.pi

        for l in range(lmax + 1):
            # Build A
            napmo_library.finite_difference_matrix(
                byref(atgrid.radial_grid), A, l)

            for m in range(-l, l + 1):
                # Build b
                b[1: -1] = sph_expansion[:, idx]
                b[1: -1] *= atgrid.radial_grid.points * pi4

                if l == 0:
                    b[0] = u_00

                # Solve
                x = np.linalg.solve(A, b)
                U[:, idx] = x[1:-1]

                idx += 1
        offset += atgrid.size

    return U

array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.finite_difference_matrix.restype = None
napmo_library.finite_difference_matrix.argtypes = [
    POINTER(RadialGrid),
    array_2d_double,
    c_int
]
