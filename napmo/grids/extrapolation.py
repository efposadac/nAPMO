# file: extrapolation.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import napmo


class Extrapolation(object):

    """
    Python interface to the Extrapolation C++ class:

    Given a set of points, can calculate the adjacent points by extrapolation.
    """

    def __init__(self):
        super(Extrapolation, self).__init__()

    def eval_left(self, x):
        return napmo.cext.Extrapolation_eval_left(self._this, x)

    def eval_right(self, x):
        return napmo.cext.Extrapolation_eval_right(self._this, x)

    def deriv_left(self, x):
        return napmo.cext.Extrapolation_deriv_left(self._this, x)

    def deriv_right(self, x):
        return napmo.cext.Extrapolation_deriv_right(self._this, x)

    def has_tail(self):
        return napmo.cext.Extrapolation_has_tail(self._this)

    def to_string(self):
        '''Return an extrapolation object in string respresentation'''
        return self.__class__.__name__


class CuspExtrapolation(Extrapolation):

    """
    Extrapolation of the form :math:`a \exp^x`
    """

    def __init__(self):
        self._this = napmo.cext.CuspExtrapolation_new()
        super(CuspExtrapolation, self).__init__()

    def prepare(self, cspline):
        napmo.cext.Extrapolation_prepare(self._this, cspline._obj)


class PowerExtrapolation(Extrapolation):

    """
    Extrapolation of the form :math:`a x^n`
    """

    def __init__(self, power):
        self._this = napmo.cext.PowerExtrapolation_new(power)
        self._power = napmo.cext.PowerExtrapolation_get_power(self._this)
        super(PowerExtrapolation, self).__init__()

    def prepare(self, cspline):
        napmo.cext.Extrapolation_prepare(self._this, cspline._obj)

    def to_string(self):
        return 'PowerExtrapolation %s' % repr(self.power)

    @property
    def power(self):
        return self._power


class PotentialExtrapolation(Extrapolation):

    """
    An extrapolation suitable for solutions of the Poisson equation.

    The prefactor of the left and right polynomial are such that the function remains
    continuous at the last grid point of the spline. The power of R on the left and right
    side is consistent with the boundary conditions of the solutions of the Poisson
    equation.

    Args:
        l (int) : The angular momentum for which the potential is computed.
    """

    def __init__(self, l):
        self.l = l
        self._this = napmo.cext.PotentialExtrapolation_new(l)
        super(PotentialExtrapolation, self).__init__()

    def prepare(self, cspline):
        napmo.cext.Extrapolation_prepare(self._this, cspline._obj)

    def to_string(self):
        return 'PotentialExtrapolation %s' % repr(self.l)
