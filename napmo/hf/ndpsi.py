# file: npsi.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


def compute_dpsi(grid, lmax, phi, doi, oi, ri, V_nuc, J, mass):
    """
    Computes \Delta \psi from Becke's paper

    Args:
        grid (BeckeGrid): Molecular grid.
        lmax (int): Limit for spherical expansion.
        phi (ndarray): Orbital function evaluated on the grid.
        doi (float): Delta orbital energy.
        oi (float): Orbital energy.
        ri (ndarray): Residual evaluated on the grid.
        V_nuc (ndarray): Nuclear potential evaluated on the grid.
        J (ndarray): Two-particles potential evaluated on the grid.
    """
    offset = 0

    dp_s = []

    _psi = (phi * doi) - ri

    _psi *= grid.becke_weights

    for atgrid in grid.atgrids:

        psi = _psi[offset:offset + atgrid.size]

        # Spherical expansion
        sph_expansion = atgrid.spherical_expansion(lmax, psi)

        V_sph = atgrid.spherical_average(
            (V_nuc[offset:offset + atgrid.size] +
             J[offset:offset + atgrid.size]) * grid.becke_weights[offset:offset + atgrid.size])

        res = []
        idx = 0
        for l in range(lmax + 1):
            for m in range(-l, l + 1):

                aux = np.array(sph_expansion[:, idx])

                # Build b
                phi = napmo.CubicSpline(aux,
                                        rtransform=atgrid.radial_grid.rtransform)

                rtf = phi.rtransform
                rgrid = napmo.RadialGrid(rtransform=rtf)
                radii = rtf.radius_all()

                fy = -2 * phi.y
                fd = -2 * phi.dx
                f = napmo.CubicSpline(fy, fd, rtf)

                b = napmo.CubicSpline(2 / radii, -2 / radii**2, rtf)

                a = napmo.CubicSpline(
                    2 * (V_sph + oi - doi + (-l * (l + 1) * radii ** -2)), rtransform=rtf)

                bcs = (0.0, None, 0.0, None)

                p = napmo.solve_ode2(
                    b, a, f, bcs, napmo.PowerExtrapolation(-l - 1))

                res.append(p)
                idx += 1

        offset += atgrid.size
        dp_s.append(res)

    dpsi = np.zeros(grid.size, dtype=np.float64)
    for i in range(grid.ncenter):
        grid.evaluate_decomposition(i, dp_s[i][:], dpsi)

    return dpsi / mass  # Check THIS!
