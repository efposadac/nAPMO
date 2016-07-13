# file: poisson_solver.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co
from __future__ import print_function

from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

import napmo


def poisson_solver(grid, dens, lmax):
    """
    Returns the spherically expanded potential :math:`U_{\ell m}(r)` obtained following Becke's procedure.

    References:
        Becke, A. D. Dickson, R. M. Numerical solution of Poisson's equation in polyatomic molecules. J. Chem. Phys. 89(5), 1988.

    Args:
        grid (BeckeGrid): Molecular grid.
        dens (ndarray): Array with the source density calculated in each point of the grid.
        lmax (int): Maximum :math:`\ell` order of the expansion.

    Returns:
        U (ndarray): Spherical expanded potential. Array with shape (nrad, lsize), where :math:`\ell_{size} = (\ell_{max} + 1)^2`
    """

    assert isinstance(grid, napmo.BeckeGrid)

    offset = 0
    U = []
    for i in range(grid.ncenter):
        atgrid = grid.atgrids[i]
        p = dens[offset:offset + atgrid.size]
        bw = grid.becke_weights[offset:offset + atgrid.size]

        # Spherical expansion
        sph_expansion = atgrid.spherical_expansion(lmax, p)

        result = []
        idx = 0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):
                aux = np.array(sph_expansion[:, idx])

                # Build b
                rho = napmo.CubicSpline(aux,
                                        rtransform=atgrid.radial_grid.rtransform)
                rtf = rho.rtransform
                rgrid = napmo.RadialGrid(rtransform=rtf)
                radii = rtf.radius_all()
                # The approach followed here is obtained after substitution of
                # u = r*V in Eq. (21) in Becke's paper. After this transformation,
                # the boundary conditions can be implemented such that the output
                # is more accurate.
                fy = -4 * np.pi * rho.y
                fd = -4 * np.pi * rho.dx
                f = napmo.CubicSpline(fy, fd, rtf)
                b = napmo.CubicSpline(2 / radii, -2 / radii**2, rtf)
                a = napmo.CubicSpline(-l * (l + 1) * radii ** -2,
                                      2 * l * (l + 1) * radii ** -3, rtf)
                # Derivation of boundary condition at rmax:
                # Multiply differential equation with r**l and integrate. Using
                # partial integration and the fact that V(r)=A/r**(l+1) for large
                # r, we find -(2l+1)A=-4pi*int_0^infty r**2 r**l rho(r) and so
                # V(rmax) = A/rmax**(l+1) = integrate(r**l
                # rho(r))/(2l+1)/rmax**(l+1)
                V_rmax = rgrid.integrate(
                    rho.y * radii**l, 1) / radii[-1]**(l + 1) / (2 * l + 1)
                # Derivation of boundary condition at rmin:
                # Same as for rmax, but multiply differential equation with r**(-l-1)
                # and assume that V(r)=B*r**l for small r.
                V_rmin = rgrid.integrate(rho.y * radii**(-l - 1),
                                         1) * radii[0]**(l) / (2 * l + 1)
                bcs = (V_rmin, None, V_rmax, None)
                v = napmo.solve_ode2(
                    b, a, f, bcs, napmo.PowerExtrapolation(-l - 1))
                result.append(v)
                idx += 1

        offset += atgrid.size
        U.append(result)

    return U


def poisson_solver_finite_differences(grid, rho, lmax):
    assert isinstance(grid, napmo.BeckeGrid)

    napmo.cext.finite_difference_matrix.restype = None
    napmo.cext.finite_difference_matrix.argtypes = [
        POINTER(napmo.RadialGrid),
        npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
        npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS'),
        npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS'),
        c_int
    ]

    offset = 0
    result = []
    for i in range(grid.ncenter):
        atgrid = grid.atgrids[i]
        p = rho[offset:offset + atgrid.size]
        bw = grid.becke_weights[offset:offset + atgrid.size]

        # Spherical expansion
        sph_expansion = atgrid.spherical_expansion(lmax, p)

        # Boundary condition
        # q = atgrid.integrate(p, bw)
        # u_00 = 0.0
        # if np.abs(q) > 1.0e-16:
        #     u_00 = np.sqrt(4.0 * np.pi * q)

        # Solve differential equation Ax = b
        U = np.zeros(sph_expansion.shape, dtype=np.float64)
        b = np.zeros(atgrid.radial_grid.size + 2, dtype=np.float64)
        data = np.zeros(atgrid.radial_grid.size * 3 + 2, dtype=np.float64)
        row = np.zeros(atgrid.radial_grid.size * 3 + 2, dtype=np.int32)
        col = np.zeros(atgrid.radial_grid.size * 3 + 2, dtype=np.int32)

        idx = 0
        pi4 = -4 * np.pi
        for l in range(lmax + 1):

            # Build A
            napmo.cext.finite_difference_matrix(
                byref(atgrid.radial_grid), data, row, col, l)
            for m in range(-l, l + 1):

                # Build b
                b[1:-1] = sph_expansion[:, idx]
                b[1:-1] *= atgrid.radial_grid.points * pi4
                if l == 0:
                    b[0] = u_00
                # Solve
                x = spsolve(
                    csc_matrix(
                        (data, (row, col)),
                        shape=(atgrid.radial_grid.size + 2,
                               atgrid.radial_grid.size + 2)),
                    b)
                U[:, idx] = x[1:-1]

                idx += 1
        result.append(U)
        offset += atgrid.size

    return result
