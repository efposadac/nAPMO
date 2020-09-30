# file: poisson_solver.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
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
                    b, a, f, bcs, napmo.PotentialExtrapolation(l))

                # v = napmo.solve_ode2(
                #     b, a, f, bcs, napmo.PowerExtrapolation(-l - 1))

                result.append(v)
                idx += 1

        offset += atgrid.size
        U.append(result)

    return U
