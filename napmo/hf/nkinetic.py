# file: nkinetic.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


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

        _psi = psi[start:end].copy()
        sph_exp = ag.spherical_expansion(
            lmax, _psi)

        idx = 0
        res = np.zeros([rtf.size, (lmax + 1)**2])
        for l in range(lmax + 1):

            num = l * (l + 1)
            pref = num * r2inv

            for m in range(-l, l + 1):

                # calculate Laplacian
                aux = np.array(sph_exp[:, idx])

                # phi'
                phi = napmo.CubicSpline(
                    aux * r, rtransform=rtf,
                    extrapolation=napmo.PowerExtrapolation(-l - 1))

                dphi = phi.deriv(r)

                # phi''
                dphi = napmo.CubicSpline(
                    dphi, rtransform=rtf,
                    extrapolation=napmo.PowerExtrapolation(-l - 1))

                d2phi = dphi.deriv(r)

                # Build equation
                res[:, idx] = ((-rinv * d2phi) + (pref * aux))
                idx += 1

        T.append(ag.evaluate_expansion(lmax, res))

    return 0.5 * np.concatenate(T)
