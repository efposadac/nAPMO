# file: cubic_spline.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import numpy.ctypeslib as npct

from napmo.system.cext import napmo_library as nl

from napmo.utilities.extrapolation import *
from napmo.grids.radial_transform import *


class CubicSpline(object):

    """
    A cubic spline object

    **Arguments:**

    y
        The function values at the 1D grid.

    **Optional arguments:**

    dx
        The derivative of the function values at the 1D grid. If not given,
        they are determined such that the second derivative of the cubic
        spline is continuous at the grid points.

    rtransform
        The transformation object that specifies the 1D grid. If not given,
        an identity transform is used

    extrapolation
        The extrapolation object that specifies the spline function outside
        the interval determined by the 1D grid. By default,
        CuspExtrapolation() is used.
    """

    def __init__(self, y, dx=None, rtransform=None, extrapolation=None):
        super(CubicSpline, self).__init__()

        self._y = y
        n = self._y.shape[0]

        # Set the rtransform
        if rtransform is None:
            self._rtransform = IndentityRadialTransform(n)
        else:
            self._rtransform = rtransform
            assert rtransform.size == n

        # use given derivatives or construct new ones.
        v = self._rtransform.deriv_all()
        if dx is None:
            self._dt = np.zeros(n, np.float64)
            nl.solve_cubic_spline_system(self._y, self._dt, n)
            self._dx = self._dt / v
        else:
            assert dx.shape[0] == n
            self._dx = dx
            self._dt = dx * v

        # Only exponential extrapolation is needed for now
        if extrapolation is None:
            self._extrapolation = CuspExtrapolation()
        else:
            self._extrapolation = extrapolation

        self._this = nl.CubicSpline_new(self._y, self._dt,
                                        self._extrapolation._this,
                                        self._rtransform._this, n)

    def __del__(self):
        nl.CubicSpline_del(self._this)

    def __call__(self, new_x, new_y=None):
        '''evaluate the spline on a grid

           **Arguments:**

           new_x
                A numpy array with the x-values at which the spline must be
                evaluated.

           **Optional arguments:**

           new_y
                When given, it is used as output argument. This array must have
                the same size of new_x.

           **Returns:** new_y
        '''
        new_n = new_x.shape[0]
        if new_y is None:
            new_y = np.zeros(new_n, np.float64)
        else:
            assert new_y.shape[0] == new_n

        nl.CubicSpline_eval(self._this, new_x, new_y, new_n)

        return new_y

    def deriv(self, new_x, new_dx=None):
        '''Evaluate the derivative of the spline (towards x) on a grid

           **Arguments:**

           new_x
                A numpy array with the x-values at which the spline must be
                evaluated.

           **Optional arguments:**

           new_dx
                When given, it is used as output argument. This array must have
                the same size of new_x.

           **Returns:** new_dx
        '''
        new_n = new_x.shape[0]
        if new_dx is None:
            new_dx = np.zeros(new_n, np.float64)
        else:
            assert new_dx.shape[0] == new_n

        nl.CubicSpline_eval_deriv(self._this, new_x, new_dx, new_n)
        return new_dx

    @property
    def y(self):
        return self._y

    @property
    def dx(self):
        return self._dx

    @property
    def first_x(self):
        return nl.CubicSpline_get_first_x(self._this)

    @property
    def last_x(self):
        return nl.CubicSpline_get_last_x(self._this)

    @property
    def extrapolation(self):
        return self._extrapolation

    @property
    def rtransform(self):
        return self._rtransform

array_1d_double = npct.ndpointer(
    dtype=np.double, ndim=1, flags='CONTIGUOUS')

nl.CubicSpline_new.restype = c_void_p
nl.CubicSpline_new.argtypes = [array_1d_double,
                               array_1d_double, c_void_p, c_void_p, c_int]

nl.CubicSpline_del.restype = None
nl.CubicSpline_del.argtypes = [c_void_p]

nl.CubicSpline_eval.restype = None
nl.CubicSpline_eval.argtypes = [
    c_void_p, array_1d_double, array_1d_double, c_int]

nl.CubicSpline_eval_deriv.restype = None
nl.CubicSpline_eval_deriv.argtypes = [
    c_void_p, array_1d_double, array_1d_double, c_int]

nl.CubicSpline_get_first_x.restype = c_double
nl.CubicSpline_get_first_x.argtypes = [c_void_p]

nl.CubicSpline_get_last_x.restype = c_double
nl.CubicSpline_get_last_x.argtypes = [c_void_p]

nl.CubicSpline_get_rtransform.restype = c_void_p
nl.CubicSpline_get_rtransform.argtypes = [c_void_p]

nl.CubicSpline_get_extrapolation.restype = c_void_p
nl.CubicSpline_get_extrapolation.argtypes = [c_void_p]

nl.solve_cubic_spline_system.restype = None
nl.solve_cubic_spline_system.argtypes = [
    array_1d_double, array_1d_double, c_int]
