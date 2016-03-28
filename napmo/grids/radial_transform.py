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

from napmo.system.cext import napmo_library as nl
from napmo.utilities.int1d import *


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

        nl.RTransform_get_npoint.restype = c_int
        nl.RTransform_get_npoint.argtypes = [c_void_p]

        self._size = nl.RTransform_get_npoint(self._this)
        self._points = np.arange(self.size, dtype=np.float64)

    def __dealloc__(self):
        if self._this != NULL:
            nl.RTransform_del.restype = None
            nl.RTransform_del.argtypes = [c_void_p]
            nl.RTransform_del(self._this)

    def radius(self, t):
        if isinstance(t, numbers.Number):
            nl.RTransform_radius.restype = c_double
            nl.RTransform_radius.argtypes = [c_void_p, c_double]
            return nl.RTransform_radius(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.radius_all(points=t)
        else:
            raise NotImplementedError

    def deriv(self, t):
        if isinstance(t, numbers.Number):
            nl.RTransform_deriv.restype = c_double
            nl.RTransform_deriv.argtypes = [c_void_p, c_double]
            return nl.RTransform_deriv(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv_all(points=t)
        else:
            raise NotImplementedError

    def deriv2(self, t):
        if isinstance(t, numbers.Number):
            nl.RTransform_deriv2.restype = c_double
            nl.RTransform_deriv2.argtypes = [c_void_p, c_double]
            return nl.RTransform_deriv2(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv2_all(points=t)
        else:
            raise NotImplementedError

    def deriv3(self, t):
        if isinstance(t, numbers.Number):
            nl.RTransform_deriv3.restype = c_double
            nl.RTransform_deriv3.argtypes = [c_void_p, c_double]
            return nl.RTransform_deriv3(self._this, t)
        elif isinstance(t, np.ndarray):
            return self.deriv3_all(points=t)
        else:
            raise NotImplementedError

    def inv(self, r):
        if isinstance(r, numbers.Number):
            nl.RTransform_inv.restype = c_double
            nl.RTransform_inv.argtypes = [c_void_p, c_double]
            return nl.RTransform_inv(self._this, r)
        elif isinstance(r, np.ndarray):
            return self.inv_all(points=r)
        else:
            raise NotImplementedError

    def radius_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        nl.RTransform_radius_array.restype = None
        nl.RTransform_radius_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        nl.RTransform_radius_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        nl.RTransform_deriv_array.restype = None
        nl.RTransform_deriv_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        nl.RTransform_deriv_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv2_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        nl.RTransform_deriv2_array.restype = None
        nl.RTransform_deriv2_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        nl.RTransform_deriv2_array(
            self._this, points, data, points.shape[0])

        return data

    def deriv3_all(self, points=None):
        if points is None:
            points = self._points

        data = np.zeros(points.shape[0], dtype=np.float64)

        nl.RTransform_deriv3_array.restype = None
        nl.RTransform_deriv3_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        nl.RTransform_deriv3_array(
            self._this, points, data, points.shape[0])

        return data

    def inv_all(self, points):
        data = np.zeros(points.shape[0], dtype=np.float64)

        nl.RTransform_inv_array.restype = None
        nl.RTransform_inv_array.argtypes = [
            c_void_p,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            c_int]

        nl.RTransform_inv_array(
            self._this, points, data, points.shape[0])

        return data

    def get_default_int1d(self):
        return SimpsonIntegrator1D()

    @property
    def size(self):
        return self._size


class IndentityRadialTransform(RadialTransform):

    """
    For testing only
    """

    def __init__(self, size):
        nl.IdentityRTransform_new.restype = c_void_p
        nl.IdentityRTransform_new.argtypes = [c_int]
        self._this = nl.IdentityRTransform_new(size)

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
        nl.PowerRTransform_new.restype = c_void_p
        nl.PowerRTransform_new.argtypes = [c_double, c_double, c_int]
        self._this = nl.PowerRTransform_new(rmin, rmax, size)

        nl.PowerRTransform_get_rmin.restype = c_double
        nl.PowerRTransform_get_rmin.argtypes = [c_void_p]
        self._rmin = nl.PowerRTransform_get_rmin(self._this)

        nl.PowerRTransform_get_rmax.restype = c_double
        nl.PowerRTransform_get_rmax.argtypes = [c_void_p]
        self._rmax = nl.PowerRTransform_get_rmax(self._this)

        nl.PowerRTransform_get_power.restype = c_double
        nl.PowerRTransform_get_power.argtypes = [c_void_p]
        self._power = nl.PowerRTransform_get_power(self._this)

        super(PowerRadialTransform, self).__init__()

    def to_string(self):
        return ' '.join(['PowerRTransform', repr(self.rmin), repr(self.rmax), repr(self.size)])

    def get_default_int1d(self):
        return StubIntegrator1D()

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
        nl.ChebyshevRTransform_new.restype = c_void_p
        nl.ChebyshevRTransform_new.argtypes = [c_double, c_int]
        self._this = nl.ChebyshevRTransform_new(radii, size)

        nl.ChebyshevRTransform_get_radii.restype = c_double
        nl.ChebyshevRTransform_get_radii.argtypes = [c_void_p]
        self._radii = nl.ChebyshevRTransform_get_radii(self._this)

        super(ChebyshevRadialTransform, self).__init__()

    def get_default_int1d(self):
        return StubIntegrator1D()

    def to_string(self):
        return ' '.join(['ChebyshevRTransform', repr(self.radii), repr(self.size)])

    @property
    def radii(self):
        return self._radii
