# file: extrapolation.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

from ctypes import *
from napmo.system.cext import napmo_library as nl


class Extrapolation(object):

    """
    docstring for Extrapolation
    """

    def __init__(self):
        super(Extrapolation, self).__init__()

    def eval_left(self, x):
        return nl.Extrapolation_eval_left(self._this, x)

    def eval_right(self, x):
        return nl.Extrapolation_eval_right(self._this, x)

    def deriv_left(self, x):
        return nl.Extrapolation_deriv_left(self._this, x)

    def deriv_right(self, x):
        return nl.Extrapolation_deriv_right(self._this, x)

    def has_tail(self):
        return nl.Extrapolation_has_tail(self._this)

    def to_string(self):
        '''Return an extrapolation object in string respresentation'''
        return self.__class__.__name__


class CuspExtrapolation(Extrapolation):

    """docstring for CuspExtrapolation"""

    def __init__(self):
        self._this = nl.CuspExtrapolation_new()
        super(CuspExtrapolation, self).__init__()

    def prepare(self, cspline):
        nl.Extrapolation_prepare(self._this, cspline._obj)


class PowerExtrapolation(Extrapolation):

    """
    docstring for PowerExtrapolation
    """

    def __init__(self, power):
        self._this = nl.PowerExtrapolation_new(power)
        self._power = nl.PowerExtrapolation_get_power(self._this)
        super(PowerExtrapolation, self).__init__()

    def prepare(self, cspline):
        nl.Extrapolation_prepare(self._this, cspline._obj)

    def to_string(self):
        return 'PowerExtrapolation %s' % repr(self.power)

    @property
    def power(self):
        return self._power

nl.Extrapolation_prepare.restype = None
nl.Extrapolation_prepare.argtypes = [c_void_p, c_void_p]

nl.Extrapolation_del.restype = None
nl.Extrapolation_del.argtypes = [c_void_p]

nl.Extrapolation_eval_left.restype = c_double
nl.Extrapolation_eval_left.argtypes = [c_void_p, c_double]

nl.Extrapolation_eval_right.restype = c_double
nl.Extrapolation_eval_right.argtypes = [c_void_p, c_double]

nl.Extrapolation_deriv_left.restype = c_double
nl.Extrapolation_deriv_left.argtypes = [c_void_p, c_double]

nl.Extrapolation_deriv_right.restype = c_double
nl.Extrapolation_deriv_right.argtypes = [c_void_p, c_double]

nl.Extrapolation_has_tail.restype = c_bool
nl.Extrapolation_has_tail.argtypes = [c_void_p]

nl.ZeroExtrapolation_new.restype = c_void_p

nl.CuspExtrapolation_new.restype = c_void_p

nl.PowerExtrapolation_new.restype = c_void_p
nl.PowerExtrapolation_new.argtypes = [c_double]

nl.PowerExtrapolation_get_power.restype = c_double
nl.PowerExtrapolation_get_power.argtypes = [c_void_p]
