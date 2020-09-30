# file: cubic_spline.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


class CubicSpline(object):

    """
    A cubic spline object

    Args:

        y : The function values at the 1D grid.
        dx (optional) : The derivative of the function values at the 1D grid. If not given,
            they are determined such that the second derivative of the cubic
            spline is continuous at the grid points.

        rtransform (optional) : The transformation object that specifies the 1D grid. If not given,
            an identity transform is used

        extrapolation (optional) : The extrapolation object that specifies the spline function outside
            the interval determined by the 1D grid. By default,
            CuspExtrapolation() is used.
    """

    def __init__(self, y, dx=None, rtransform=None, extrapolation=None):
        super(CubicSpline, self).__init__()

        self._y = y
        n = self._y.shape[0]

        # Set the rtransform
        if rtransform is None:
            self._rtransform = napmo.IndentityRadialTransform(n)
        else:
            self._rtransform = rtransform
            assert rtransform.size == n

        # use given derivatives or construct new ones.
        v = self._rtransform.deriv_all()
        if dx is None:
            self._dt = np.zeros(n, np.float64)
            napmo.cext.solve_cubic_spline_system(self._y, self._dt, n)
            self._dx = self._dt / v
        else:
            assert dx.shape[0] == n
            self._dx = dx
            self._dt = dx * v

        # Only exponential extrapolation is needed for now
        if extrapolation is None:
            self._extrapolation = napmo.CuspExtrapolation()
        else:
            self._extrapolation = extrapolation

        self._this = napmo.cext.CubicSpline_new(self._y, self._dt,
                                                self._extrapolation._this,
                                                self._rtransform._this, n)

    def __del__(self):
        napmo.cext.CubicSpline_del(self._this)

    def __call__(self, new_x, new_y=None):
        """
        evaluate the spline on a grid

        Args:
            new_x (ndarray): A numpy array with the x-values at which the spline must be evaluated.
            new_y (ndarray, optional) : When given, it is used as output argument. This array must have
                the same size of new_x.
        Returns:
            new_y
        """
        new_n = new_x.shape[0]
        if new_y is None:
            new_y = np.zeros(new_n, np.float64)
        else:
            assert new_y.shape[0] == new_n

        napmo.cext.CubicSpline_eval(self._this, new_x, new_y, new_n)

        return new_y

    def deriv(self, new_x, new_dx=None):
        """
        Evaluate the derivative of the spline (towards x) on a grid

        Args:

            new_x (ndarray) : A numpy array with the x-values at which the spline must be
                evaluated.
            new_dx (ndarray, optional) : When given, it is used as output argument. This array must have
                the same size of new_x.
        Returns:
            new_dx
        """
        new_n = new_x.shape[0]
        if new_dx is None:
            new_dx = np.zeros(new_n, np.float64)
        else:
            assert new_dx.shape[0] == new_n

        napmo.cext.CubicSpline_eval_deriv(self._this, new_x, new_dx, new_n)
        return new_dx

    @property
    def y(self):
        return self._y

    @property
    def dx(self):
        return self._dx

    @property
    def first_x(self):
        return napmo.cext.CubicSpline_get_first_x(self._this)

    @property
    def last_x(self):
        return napmo.cext.CubicSpline_get_last_x(self._this)

    @property
    def extrapolation(self):
        return self._extrapolation

    @property
    def rtransform(self):
        return self._rtransform
