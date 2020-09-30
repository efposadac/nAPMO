# file: npsi.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


def compute_dpsi(grid, lmax, phi, doi, oi, ri, V_tot, mass_inv, charge):
    """
    Computes \Delta \psi from Becke's paper

    Args:
        grid (BeckeGrid): Molecular grid.
        lmax (int): Limit for spherical expansion.
        phi (ndarray): Orbital function evaluated on the grid.
        doi (float): Delta orbital energy.
        oi (float): Orbital energy.
        ri (ndarray): Residual evaluated on the grid.
        V_tot (ndarray): Electrostatic potential evaluated on the grid.
        mass_inv (float) : Mass inverse
    """
    offset = 0

    dp_s = []

    _psi = (phi * doi) - ri

    _psi *= grid.becke_weights

    for atgrid in grid.atgrids:

        psi = _psi[offset:offset + atgrid.size].copy()

        # Spherical expansion
        sph_expansion = atgrid.spherical_expansion(lmax, psi)

        V_sph = atgrid.spherical_average(
            V_tot[offset:offset + atgrid.size] * grid.becke_weights[offset:offset + atgrid.size])

        rtf = atgrid.radial_grid.rtransform

        r = rtf.radius_all()
        r2 = r * r
        rinv = 1.0 / r
        r2inv = 1.0 / r2

        res = []
        idx = 0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):

                aux = np.array(sph_expansion[:, idx])

                # Build b
                mu = napmo.CubicSpline(aux,
                                       rtransform=atgrid.radial_grid.rtransform)

                fy = -2.0 / mass_inv * mu.y
                fd = -2.0 / mass_inv * mu.dx

                f = napmo.CubicSpline(fy, fd, rtf)

                b = napmo.CubicSpline(2.0 * rinv, -2.0 * r2inv, rtf)

                a = napmo.CubicSpline(-2.0 / mass_inv * ((l * (l + 1.0) * r2inv) + (V_sph - oi - doi)), rtransform=rtf)

                bcs = (0.0, None, 0.0, None)

                p = napmo.solve_ode2(
                    b, a, f, bcs, napmo.PotentialExtrapolation(l))

                res.append(p)
                idx += 1

        offset += atgrid.size
        dp_s.append(res)

    dpsi = np.zeros(grid.size, dtype=np.float64)
    for i in range(grid.ncenter):
        grid.evaluate_decomposition(i, dp_s[i][:], dpsi)

    return dpsi
