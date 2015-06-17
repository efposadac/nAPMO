# file: analytical_integration.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
import numpy as np

from math import *


def obaraSaika_recursion(PA, PB, gamma, l_a, l_b):
    """
    Perform the Obara-Saika (1988) recursion to calculate the overlap integral:

    :math:`<\phi_A|\phi_B>`

    Where each :math:`\phi` corresponds to a GTO PrimitiveGaussian

    Args:
        PA (numpy.ndarray(3)) : Origin of :math:`\phi_A`
        PB (numpy.ndarray(3)) : Origin of :math:`\phi_B`
        gamma (float64) : Reduced exponent (see: Gaussian product theorem)
        l_a (int) : Angular moment index of :math:`\phi_A`
        l_b (int) : Angular moment index of :math:`\phi_B`

    Returns:
        x, y, z (tuple) : x, y, and z components of the recursion.
    """
    pp = 1./(2.*gamma)

    max_l = max(l_a, l_b)

    x = np.zeros([max_l+2, max_l+2])
    y = np.zeros([max_l+2, max_l+2])
    z = np.zeros([max_l+2, max_l+2])

    x[0, 0] = 1.0
    y[0, 0] = 1.0
    z[0, 0] = 1.0

    # Upward recursion in j for i=0
    x[0, 1] = PB[0]
    y[0, 1] = PB[1]
    z[0, 1] = PB[2]

    for j in range(1, l_b):
        x[0, j+1] = PB[0]*x[0, j]
        y[0, j+1] = PB[1]*y[0, j]
        z[0, j+1] = PB[2]*z[0, j]
        x[0, j+1] = x[0, j+1] + j*pp*x[0, j-1]
        y[0, j+1] = y[0, j+1] + j*pp*y[0, j-1]
        z[0, j+1] = z[0, j+1] + j*pp*z[0, j-1]

    # Upward recursion in i for all j
    x[1, 0] = PA[0]
    y[1, 0] = PA[1]
    z[1, 0] = PA[2]

    for j in range(1, l_b + 1):
        x[1, j] = PA[0]*x[0, j]
        y[1, j] = PA[1]*y[0, j]
        z[1, j] = PA[2]*z[0, j]
        x[1, j] = x[1, j] + j*pp*x[0, j-1]
        y[1, j] = y[1, j] + j*pp*y[0, j-1]
        z[1, j] = z[1, j] + j*pp*z[0, j-1]

    for i in range(1, l_a):
        x[i+1, 0] = PA[0]*x[i, 0]
        y[i+1, 0] = PA[1]*y[i, 0]
        z[i+1, 0] = PA[2]*z[i, 0]
        x[i+1, 0] = x[i+1, 0] + i*pp*x[i-1, 0]
        y[i+1, 0] = y[i+1, 0] + i*pp*y[i-1, 0]
        z[i+1, 0] = z[i+1, 0] + i*pp*z[i-1, 0]
        for j in range(1, l_b + 1):
            x[i+1, j] = PA[0]*x[i, j]
            y[i+1, j] = PA[1]*y[i, j]
            z[i+1, j] = PA[2]*z[i, j]
            x[i+1, j] = x[i+1, j] + i*pp*x[i-1, j]
            y[i+1, j] = y[i+1, j] + i*pp*y[i-1, j]
            z[i+1, j] = z[i+1, j] + i*pp*z[i-1, j]
            x[i+1, j] = x[i+1, j] + j*pp*x[i, j-1]
            y[i+1, j] = y[i+1, j] + j*pp*y[i, j-1]
            z[i+1, j] = z[i+1, j] + j*pp*z[i, j-1]

    return x, y, z


def kronecker_delta(a, b):
    """
    Calculates the delta of Kronecker for `a` and `b`.

    :math:`\delta_{ij} = 0`  if :math:`i \\neq j`, otherwise, :math:`\delta_{ij} = 1`

    Args:
        a (int) : i in  :math:`\delta_{ij}`
        b (int) : j in  :math:`\delta_{ij}`
    """
    aa = np.abs(a)
    bb = np.abs(b)
    output = int((float((aa+bb+2)-np.abs(aa-bb))) / (float((aa+bb+2) + np.abs(aa-bb))))
    return output
