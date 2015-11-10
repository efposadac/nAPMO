# file: chebyshev.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np


def chebgauss(n):
    """
    Computes abscissas and weights for the Gauss-Chebyshev quadrature of second kind.
    """
    r = np.zeros([n], dtype=np.float64)
    w = np.zeros([n], dtype=np.float64)

    for i in range(n):
        r[i] = np.cos(np.pi * (i + 1) / (n + 1.0))
        w[i] = (np.pi / (n + 1.0) * (1.0 - r[i]**2))
    return r, w


def chebgauss_transformed(n):
    """Computes the abscissas and weights for transformed Gauss-Chebyshev quadrature of second kind.
    Special case for the integral :math:`\int f(x) dx`. since Gauss-Chebyshev is for :math:`\int f(x) \sqrt{1-x^2}`.

    References:
        Perez-Jorda, J., San-Fabian, E. & Moscardo, F. A simple, reliable and efficient scheme for \
automatic numerical integration. Comput. Phys. Commun. 70, 271-284 (1992).

    Args:
        n (integer): Order of the quadrature > 0.

    Returns:
        rq (numpy.ndarray[n, 2]): Array with the abscissas and weights. r = [:, 0], w = [:, 1].
    """

    assert n > 0

    rq = np.zeros([n, 2], order='F', dtype=np.float64)

    count = 1
    for i in range(int((n + 1) / 2)):
        d = n + 1.
        s = np.sin((count * np.pi) / d)
        c = np.cos((count * np.pi) / d)
        s2 = s * s
        s4 = s2 * s2
        rq[i, 0] = 1. + (2. / np.pi) * (
            1. + (2. / 3.) * s2) * c * s - (2. * count) / d
        rq[n - i - 1, 0] = -rq[i, 0]
        rq[i, 1] = 16. / (3. * d) * s4
        rq[n - i - 1, 1] = rq[i, 1]
        count += 1

    return rq


def chebgauss_transformed_integrate(f, eps=1.0e-10, max_iter=10):
    """Computes the integral of f(x) by using transformed Gauss-Chebyshev quadrature of second kind.

    References:

        Perez-Jorda, J., San-Fabian, E. & Moscardo, F. A simple, reliable and efficient scheme for \
automatic numerical integration. Comput. Phys. Commun. 70, 271-284 (1992).

    Args:
        f (function, auto-callable): Function :math:`f(x)`, to be integrate.
        eps(float64, optional): Tolerance, Default is :math:`10^{-10}`
        max_iter (int, optional): Maximum number of iterations. Default is :math:`10`

    Returns:
        tuple: Integral value and error (out, err).
    """

    n = 5
    err = np.float64(1.0)
    count = 0
    aux = np.float64(0.0)

    while err > eps:
        rq = chebgauss_transformed(n)
        out = np.float64(0.0)
        for i in range(n):
            out += f(rq[i, 0]) * rq[i, 1]

        err = np.abs(out - aux)
        count += 1
        if count > max_iter:
            print('Convergence failed!')
            break
        n += n + 1
        aux = out

    return out, err