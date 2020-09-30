# file: nkinetic.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo
from scipy.interpolate import splrep, splev


def compute_kinetic(grid, psi, lmax):
    """
    Computes the action of the Kinetic operator on a trial function

    :math:`-\dfrac{1}{2} \nabla^{2} \phi_{i}`

    Args:
        grid (BeckeGrid): Molecular grid
        psi (ndarray): Orbital calculated on the grid
        lmax (int): Size of the spherical expansion

    Return:
        T (ndaray): Action of :math:`-\dfrac{1}{2} \nabla^2` over ``psi``
    """
    T = []
    for i, ag in enumerate(grid.atgrids):
        start = i * ag.size
        end = start + ag.size

        # Calculate spherical expansion
        rtf = ag.radial_grid.rtransform

        r = rtf.radius_all()
        r2 = r * r
        rinv = 1.0 / r
        r2inv = 1.0 / r2

        _psi = psi[start:end].copy() * grid.becke_weights[start:end]
        sph_exp = ag.spherical_expansion(lmax, _psi)

        idx = 0
        res = []

        for l in range(lmax + 1):

            num = l * (l + 1)
            pref = num * r2inv

            for m in range(-l, l + 1):

                # calculate Laplacian
                aux = np.array(sph_exp[:, idx]).copy()

                phi = splrep(r, aux * r, k=5)
                d2phi = splev(r, phi, der=2)

                ####################################################
                # phi'
                # phi = napmo.CubicSpline(
                #     aux * r, rtransform=rtf,
                #     extrapolation=napmo.PowerExtrapolation(-l - 1))

                # dphi = phi.deriv(r)

                # # phi''
                # dphi = napmo.CubicSpline(
                #     dphi, rtransform=rtf,
                #     extrapolation=napmo.PowerExtrapolation(-l - 1))

                # d2phi = dphi.deriv(r)
                ####################################################

                # Build equation
                res.append(napmo.CubicSpline(
                    ((rinv * d2phi) - (pref * aux)), rtransform=rtf,
                    extrapolation=napmo.PowerExtrapolation(-l - 1)))

                idx += 1

        T.append(res)

    with napmo.runtime.timeblock('Interpolation'):
        U = np.zeros(grid.size, dtype=np.float64)
        for i in range(grid.ncenter):
            grid.evaluate_decomposition(i, T[i][:], U)

    return -0.5 * U
