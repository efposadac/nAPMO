# file : radial_transform.py
# nAPMO package
# Copyright(c)2014, Edwin Fernando Posada
# All rights reserved.
# Version : 0.1
# efposadac @sissa.it

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import numpy.ctypeslib as npct
import numbers

import napmo


class RadialTransform(object):

    """
    A definition of (radial) grid points by means of a transformation.

       The definition starts from a uniform 1D grid with spacing 1 and starting
       point 0: 0, 1, 2, 3, ... npoint-1. These values are defined on the
       so-called t-axis. The transformation is a function r=f(t) that defines
       the actual grid points on the r-axis: f(0), f(1), f(2), ... f(npoint-1).
       Different implementation for the function f are available.

    Abstract class for radial transformation. Not use it directly
    """

    def __init__(self):
        super(RadialTransform, self).__init__()

        napmo.cext.RTransform_get_npoint.restype = c_int
        napmo.cext.RTransform_get_npoint.argtypes = [c_void_p]

        self._size = napmo.cext.RTransform_get_npoint(self._this)
        self._points = np.arange(self.size, dtype=np.float64)

    def __dealloc__(self):
        if self._this != NULL:
            napmo.cext.RTransform_del.restype = None
            napmo.cext.RTransform_del.argtypes = [c_void_p]
            napmo.cext.RTransform_del(self._this)

    def radius(self, t):
        if isinstance(t, numbers.Number):
            napmo.cext.RTransform_radius.restype = c_double
            napmo.cext.RTransform_radius.argtypes = [c_void_p, c_double]
            return napmo.cext.RTransform_radius(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.radius_all(points=t)
        else:
            raise NotImplementedError

    def deriv(self, t):
        if isinstance(t, numbers.Number):
            napmo.cext.RTransform_deriv.restype = c_double
            napmo.cext.RTransform_deriv.argtypes = [c_void_p, c_double]
            return napmo.cext.RTransform_deriv(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv_all(points=t)
        else:
            raise NotImplementedError

    def deriv2(self, t):
        if isinstance(t, numbers.Number):
            napmo.cext.RTransform_deriv2.restype = c_double
            napmo.cext.RTransform_deriv2.argtypes = [c_void_p, c_double]
            return napmo.cext.RTransform_deriv2(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv2_all(points=t)
        else:
            raise NotImplementedError

    def deriv3(self, t):
        if isinstance(t, numbers.Number):
            napmo.cext.RTransform_deriv3.restype = c_double
            napmo.cext.RTransform_deriv3.argtypes = [c_void_p, c_double]
            return napmo.cext.RTransform_deriv3(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv3_all(points=t)
        else:
            raise NotImplementedError

    def inv(self, r):
        if isinstance(r, numbers.Number):
            napmo.cext.RTransform_inv.restype = c_double
            napmo.cext.RTransform_inv.argtypes = [c_void_p, c_double]
            return napmo.cext.RTransform_inv(self._this, r)
        elif isinstance(r, np.ndarray):
            return self.inv_all(points=r)
        else:
            raise NotImplementedError

    def radius_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        napmo.cext.RTransform_radius_array.restype = None
        napmo.cext.RTransform_radius_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        napmo.cext.RTransform_radius_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        napmo.cext.RTransform_deriv_array.restype = None
        napmo.cext.RTransform_deriv_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        napmo.cext.RTransform_deriv_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv2_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        napmo.cext.RTransform_deriv2_array.restype = None
        napmo.cext.RTransform_deriv2_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        napmo.cext.RTransform_deriv2_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv3_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        napmo.cext.RTransform_deriv3_array.restype = None
        napmo.cext.RTransform_deriv3_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        napmo.cext.RTransform_deriv3_array(
            self._this, points, data, points.shape[0])

        return data

    def inv_all(self, points):
        data = np.zeros(points.shape[0], dtype=np.float64)

        napmo.cext.RTransform_inv_array.restype = None
        napmo.cext.RTransform_inv_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        napmo.cext.RTransform_inv_array(
            self._this, points, data, points.shape[0])

        return data

    def get_default_int1d(self):
        return napmo.SimpsonIntegrator1D()

    @property
    def size(self):
        return self._size


class IndentityRadialTransform(RadialTransform):

    """
    For testing only
    """

    def __init__(self, size):
        napmo.cext.IdentityRTransform_new.restype = c_void_p
        napmo.cext.IdentityRTransform_new.argtypes = [c_int]
        self._this = napmo.cext.IdentityRTransform_new(size)

        super(IndentityRadialTransform, self).__init__()

    def to_string(self):
        return ' '.join(['IdentityRTransform', repr(self.size)])


class PowerRadialTransform(RadialTransform):

    """
    A power grid.

       The grid points are distributed as follows:

       .. math:: r_i = r_0 i^{\alpha}

       with

       .. math::
            \alpha = \frac{\ln r_{N-1} - \ln r_0}{\ln N-1}
    """

    def __init__(self, rmin, rmax, size):
        napmo.cext.PowerRTransform_new.restype = c_void_p
        napmo.cext.PowerRTransform_new.argtypes = [c_double, c_double, c_int]
        self._this = napmo.cext.PowerRTransform_new(rmin, rmax, size)

        napmo.cext.PowerRTransform_get_rmin.restype = c_double
        napmo.cext.PowerRTransform_get_rmin.argtypes = [c_void_p]
        self._rmin = napmo.cext.PowerRTransform_get_rmin(self._this)

        napmo.cext.PowerRTransform_get_rmax.restype = c_double
        napmo.cext.PowerRTransform_get_rmax.argtypes = [c_void_p]
        self._rmax = napmo.cext.PowerRTransform_get_rmax(self._this)

        napmo.cext.PowerRTransform_get_power.restype = c_double
        napmo.cext.PowerRTransform_get_power.argtypes = [c_void_p]
        self._power = napmo.cext.PowerRTransform_get_power(self._this)

        super(PowerRadialTransform, self).__init__()

    def to_string(self):
        return ' '.join(['PowerRTransform', repr(self.rmin), repr(self.rmax), repr(self.size)])

    def get_default_int1d(self):
        return napmo.StubIntegrator1D()

    @property
    def rmin(self):
        return self._rmin

    @property
    def rmax(self):
        return self._rmax

    @property
    def power(self):
        return self._power


class ChebyshevRadialTransform(RadialTransform):

    """
    ...... something
    """

    def __init__(self, radii, size):
        napmo.cext.ChebyshevRTransform_new.restype = c_void_p
        napmo.cext.ChebyshevRTransform_new.argtypes = [c_double, c_int]
        self._this = napmo.cext.ChebyshevRTransform_new(radii, size)

        napmo.cext.ChebyshevRTransform_get_radii.restype = c_double
        napmo.cext.ChebyshevRTransform_get_radii.argtypes = [c_void_p]
        self._radii = napmo.cext.ChebyshevRTransform_get_radii(self._this)

        super(ChebyshevRadialTransform, self).__init__()

    def get_default_int1d(self):
        return napmo.StubIntegrator1D()

    def to_string(self):
        return ' '.join(['ChebyshevRTransform', repr(self.radii), repr(self.size)])

    @property
    def radii(self):
        return self._radii
