# file: angular_quadratures.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np
import numpy.ctypeslib as npct

from ctypes import *
from napmo.interfaces.c_binding import napmo_library


def lebedev(n, spherical=True):
    """Computes the Lebedev points and weights for spherical integration.

    References:
        V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59, No. 3, 477 (1999)

    Args:
        n (int): Number of angular points.

    Returns:
        tupple (numpy.ndarray): Array with the coordinates and weights (theta[n], phi[n], w[n]).

    Raises:
        ValueError: If ``n`` is not supported.
    """

    assert n > 0

    lebpoints = [
        6, 14, 110, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454,
        1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810
    ]

    if n not in lebpoints:
        print("Info: Number of angular points supported:...")
        print(lebpoints)
        raise ValueError('Invalid number of angular points', n)

    if spherical:
        t = np.empty(n)
        p = np.empty(n)
        w = np.empty(n)
        napmo_library.lebedev_spherical(n, t, p, w)
        return t, p, w
    else:
        w = np.empty(n)
        p = np.empty([n, 3])
        napmo_library.lebedev_cartesian(n, p, w)
        return p, w


def lebedev_get_order(n):
    """
    Return the order of the Lebedev quadrature from the number of points ``n``.
    """
    lebpoints = {
        6: 3, 14: 5, 110: 17, 146: 19, 170: 21, 194: 23, 230: 25, 266: 27, 302: 29, 350: 31, 434: 35,
        590: 41, 770: 47, 974: 53, 1202: 59, 1454: 65, 1730: 71, 2030: 77, 2354: 83,
        2702: 89, 3074: 95, 3470: 101, 3890: 107, 4334: 113, 4802: 119, 5294: 125, 5810: 131
    }

    return lebpoints[n]


def lebedev_integrate(func, n, args=()):
    """
    Integrate the ``func`` over unit sphere :math:`4\pi` using Lebedev quadrature of ``n`` points.

    Args:
        func (callable): A Python function or method of at least two variables: theta and phi
        n (int): Number of Lebedev points to be used in the calculation.
        args (sequence, optional): Extra arguments to pass to func.
    """
    theta, phi, w = lebedev(n)
    integral = (func(theta, phi, *args) * w).sum()
    return integral * 4.0 * np.pi


array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.lebedev_spherical.restype = None
napmo_library.lebedev_spherical.argtypes = [
    c_int,
    array_1d_double,
    array_1d_double,
    array_1d_double
]

napmo_library.lebedev_cartesian.restype = None
napmo_library.lebedev_cartesian.argtypes = [
    c_int,
    array_2d_double,
    array_1d_double
]
