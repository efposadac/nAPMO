# -*- coding: utf-8 -*-
# file: ode2.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2015 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#

'''
Finite-element second-order ODE solver
'''


from scipy.sparse.linalg import spsolve
from scipy.sparse import csc_matrix

import numpy as np
from ctypes import *

import napmo


def solve_ode2(b, a, f, bcs, extrapolation=None):
    """
    Solve a second order ODE.

    Args:

        b, a, f (CubicSpline) : Cubic splines for the given functions in the second order ODE. (See
            build_neumann for details.) These cubic splines must have identical
            RTransform objects.

        bcs : The boundary conditions (See build_neumann for details.)

        extrapolation (Extrapolation, optional) : The extrapolation object for the returned cubic spline.

    Returns:
        output (CubicSpline) : a cubic spline object with the solution that uses the same
            RTransform object as the input functions a, b and f.
    """

    def merge(y, d):
        '''Put y and d in one vector'''
        return np.array([y, d]).T.ravel()

    # Parse args.
    rtf = b.rtransform
    if rtf.to_string() != a.rtransform.to_string():
        raise ValueError('The RTransform objects of b and a do not match.')
    if rtf.to_string() != f.rtransform.to_string():
        raise ValueError('The RTransform objects of b and f do not match.')

    # Transform the given functions to the linear coordinate.
    j1 = rtf.deriv_all()
    j2 = rtf.deriv2_all()
    j3 = rtf.deriv3_all()
    j1sq = j1 * j1

    by_new = j1 * b.y - j2 / j1
    bd_new = j2 * b.y + j1sq * b.dx + (j2 * j2 - j1 * j3) / j1sq

    ay_new = a.y * j1sq
    ad_new = (a.dx * j1sq + 2 * a.y * j2) * j1

    fy_new = f.y * j1sq
    fd_new = (f.dx * j1sq + 2 * f.y * j2) * j1

    # Transform the boundary conditions
    new_bcs = (
        bcs[0], None if bcs[1] is None else bcs[1] * j1[0],
        bcs[2], None if bcs[3] is None else bcs[3] * j1[-1],
    )

    # Prepare data
    npoint = by_new.shape[0]
    b = merge(by_new, bd_new)
    a = merge(ay_new, ad_new)
    f = merge(fy_new, fd_new)

    # parse boundary conditions argument
    c_bcs = np.array([np.NaN if bcs is None else bcs
                      for bcs in new_bcs], dtype=np.float64)

    # prepare output
    coeffs = np.zeros((2 * npoint, 2 * npoint), dtype=np.float64)
    rhs = np.zeros(2 * npoint, dtype=np.float64)

    # call c routine
    with napmo.runtime.timeblock('Prepare ODE'):
        napmo.cext.build_ode2(b, a, f, c_bcs, coeffs, rhs, npoint)

    with napmo.runtime.timeblock('Solve ODE'):
        solution = spsolve(csc_matrix(coeffs), rhs)

    uy_new = solution[::2]
    ud_new = solution[1::2]

    # Transform solution back to the original coordinate.
    uy_orig = uy_new.copy()  # A copy of is needed to obtain contiguous arrays.
    ud_orig = ud_new / j1

    return napmo.CubicSpline(uy_orig, ud_orig, rtf, extrapolation)
