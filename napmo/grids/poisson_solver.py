# file: poisson_solver.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import print_function

from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

import numpy as np
from ctypes import *
import napmo


def poisson_solver(grid, _dens, lmax, sph_exp=None):
    """
    Returns the spherically expanded potential :math:`U_{\ell m}(r)` obtained following Becke's procedure.
    (finite elements).

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

    dens = _dens * grid.becke_weights
    offset = 0
    PI_4 = 4.0 * np.pi

    U = []
    for i in range(grid.ncenter):
        atgrid = grid.atgrids[i]
        p = dens[offset:offset + atgrid.size].copy()

        # Spherical expansion
        if not sph_exp:
            sph_expansion = atgrid.spherical_expansion(lmax, p)
        else:
            sph_expansion = sph_exp[i]

        result = []
        idx = 0

        rgrid = atgrid.radial_grid
        radii = rgrid.points
        rtf = rgrid.rtransform

        Dradii = rtf.deriv_all()
        D2radii = rtf.deriv2_all()
        D3radii = rtf.deriv3_all()

        radii2 = radii * radii
        radii3 = radii2 * radii

        b = napmo.CubicSpline(2 / radii, -2 / radii2, rtransform=rtf, v=Dradii)

        for l in range(lmax + 1):
            radiiL = radii**l
            radiiLm1 = radii**(-l - 1)
            l2l = l * (l + 1)
            radiiVmax = 1.0 / radii[-1]**(l + 1) / (2 * l + 1)
            radiiVmin = radii[0]**(l) / (2 * l + 1)

            # Build a
            a = napmo.CubicSpline(-l2l / radii2, 2 * l2l / radii3, rtransform=rtf, v=Dradii)

            for m in range(-l, l + 1):
                aux = np.array(sph_expansion[:, idx])

                # Build b
                rho = napmo.CubicSpline(aux, rtransform=rtf, v=Dradii)

                # The approach followed here is obtained after substitution of
                # u = r*V in Eq. (21) in Becke's paper. After this transformation,
                # the boundary conditions can be implemented such that the output
                # is more accurate.
                fy = -PI_4 * rho.y
                fd = -PI_4 * rho.dx
                f = napmo.CubicSpline(fy, fd, rtransform=rtf, v=Dradii)

                # Derivation of boundary condition at rmax:
                # Multiply differential equation with r**l and integrate. Using
                # partial integration and the fact that V(r)=A/r**(l+1) for large
                # r, we find -(2l+1)A=-4pi*int_0^infty r**2 r**l rho(r) and so
                # V(rmax) = A/rmax**(l+1) = integrate(r**l
                # rho(r))/(2l+1)/rmax**(l+1)
                V_rmax = rgrid.integrate(rho.y * radiiL, 1) * radiiVmax

                # Derivation of boundary condition at rmin:
                # Same as for rmax, but multiply differential equation with r**(-l-1)
                # and assume that V(r)=B*r**l for small r.
                V_rmin = rgrid.integrate(rho.y * radiiLm1, 1) * radiiVmin

                bcs = (V_rmin, None, V_rmax, None)

                # Solve the differential equation
                v = napmo.solve_ode2(
                    b, a, f, bcs, napmo.PotentialExtrapolation(l), j1=Dradii, j2=D2radii, j3=D3radii
                )

                result.append(v)
                idx += 1

        offset += atgrid.size
        U.append(result)

    return U
