# file: nkinetic.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


def compute_kinetic(grid, psi, lmax):
    """
    Computes the action of the Kinetic operator on a trial function

    :math:`\dfrac{1}{2} \\nabla^{2} \phi_{i}`

    Args:
        grid (BeckeGrid): Molecular grid
        psi (ndarray): Orbital calculated on the grid
        lmax (int): Size of the spherical expansion

    Return:
        T (ndaray): Action of :math:`-\dfrac{1}{2} \\nabla^2` over ``psi``
    """
    T = []
    for i, ag in enumerate(grid.atgrids):
        start = i * ag.size
        end = start + ag.size

        # Calculate spherical expansion
        rtf = ag.radial_grid.rtransform

        r = ag.radial_grid.points
        r2 = r * r
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

                d2phi = ag.radial_grid.deriv2(aux)

                # Build equation
                res.append(napmo.CubicSpline(
                    (d2phi - (pref * aux)), rtransform=rtf,
                    extrapolation=napmo.PotentialExtrapolation(l)))

                idx += 1

        T.append(res)

    with napmo.runtime.timeblock('Interpolation'):
        U = np.zeros(grid.size, dtype=np.float64)
        for i in range(grid.ncenter):
            grid.evaluate_decomposition(i, T[i][:], U)

    return -0.5 * U
