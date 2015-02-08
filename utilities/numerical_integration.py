# file: numerical_integration.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it
from __future__ import division
import numpy as np
import utilities.lebedev as lebedev


def chebgauss_t_q_w(n, rescale=False):
    """Computes the abscissas and weights for transformed Gauss-Chebyshev quadrature of second kind.

    Special case for the integral I f(x) dx. since Gauss-Chebyshev is for I f(x) sqrt(1-x^2).

    see: Perez-Jorda, J., San-Fabian, E. & Moscardo, F. A simple, reliable and efficient scheme for
    automatic numerical integration. Comput. Phys. Commun. 70, 271-284 (1992).

    returns the weights and abscissas to perform integration
    """

    assert n > 0
    assert isinstance(rescale, bool)

    w = np.zeros(n, dtype=np.float64)
    r = np.zeros(n, dtype=np.float64)

    count = 1
    for i in xrange(int((n+1)/2)):
        d = n + 1.
        s = np.sin((count*np.pi)/d)
        c = np.cos((count*np.pi)/d)
        s2 = s*s
        s4 = s2 * s2
        r[i] = 1. + (2./np.pi) * (1. + (2./3.) * s2) * c * s - (2. * count)/d
        r[n-i-1] = -r[i]
        w[i] = 16./(3.*d) * s4
        w[n-i-1] = w[i]
        count += 1

    if rescale:
        w = (w/(1.-r))/np.log(2.)
        r = 1. - np.log(1.-r)/np.log(2.)

    return r, w


def chebgauss_integrate(f, eps=1.0e-10, max_iter=10):
    """Computes the integral of f(x) by using Gauss-Chebyshev quadrature of second kind.

    see: see: Perez-Jorda, J., San-Fabian, E. & Moscardo, F. A simple, reliable and efficient scheme for
    automatic numerical integration. Comput. Phys. Commun. 70, 271-284 (1992).

    f: function f(x), to be integrate.
    eps: tolerance
    max_iter: Maximum number of iterations
    return: output type tuple with integral value and error.
    """

    n = 5
    err = np.float64(1.0)
    count = 0
    aux = np.float64(0.0)

    while err > eps:
        r, w = chebgauss_t_q_w(n)
        out = np.float64(0.0)
        for i in xrange(n):
            out += f(r[i]) * w[i]

        err = np.abs(out - aux)
        count += 1
        if count > max_iter:
            print 'Convergence failed!'
            break
        n = n + n + 1
        aux = out

    return out, err


def lebedev_q_w(n):
    """Computes the Lebedev points and weights for spherical integration.

    see: V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59, No. 3, 477 (1999)

    x, y, z: Cartesian coordinates
    w: weights
    """

    assert n > 0

    x = np.zeros(n, dtype=np.float64)
    y = np.zeros(n, dtype=np.float64)
    z = np.zeros(n, dtype=np.float64)
    w = np.zeros(n, dtype=np.float64)

    if lebedev.lebedev_compute(x, y, z, w, n) != 0:
        return x, y, z, w
    else:
        print "Info: Number of angular points supported:..."
        print "Info: 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770"
        raise ValueError('Invalid number of angular points', n)
