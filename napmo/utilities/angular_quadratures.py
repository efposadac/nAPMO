# file: angular_quadratures.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
import numpy as np
from ctypes import *
from napmo.interfaces.c_binding import *


def lebedev(n):
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

    size = (c_double * n)()
    _t = cast(size, POINTER(c_double))

    size = (c_double * n)()
    _p = cast(size, POINTER(c_double))

    size = (c_double * n)()
    _w = cast(size, POINTER(c_double))

    napmo_library.lebedev(n, _t, _p, _w)

    t = np.zeros([n], dtype=np.float64)
    p = np.zeros([n], dtype=np.float64)
    w = np.zeros([n], dtype=np.float64)

    for i in range(n):
        t[i] = _t[i]
        p[i] = _p[i]
        w[i] = _w[i]

    return t, p, w


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
