# file: radial_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

import napmo


class RadialGrid(Structure):

    '''
    An integration grid for the radial component of a spherical coordinate system
    '''

    _fields_ = [
        ("_size", c_int), ("_radii", c_double), ("_points", POINTER(c_double)),
        ("_weights", POINTER(c_double))
    ]

    def __init__(self, size=None, atomic_symbol=None, rtransform=None, int1d=None):

        if rtransform is None:
            self._size = size
            self._radii = napmo.AtomicElementsDatabase()[atomic_symbol][
                'atomic_radii_2']

            # tmp = napmo.RadialGridCheb(size, atomic_symbol)

            # self._rtransform = PowerRadialTransform(tmp.points.min(),
            #                                         tmp.points.max(),
            #                                         tmp.size)
            self._rtransform = napmo.ChebyshevRadialTransform(
                self._radii, size)
        else:
            self._rtransform = rtransform
            self._size = self._rtransform.size

        if int1d is None:
            self._int1d = self._rtransform.get_default_int1d()
        else:
            self._int1d = int1d

        self.points = self._rtransform.radius_all()
        self._points = np.ctypeslib.as_ctypes(self.points)

        self.weights = (4 * np.pi) * (
            self._rtransform.deriv_all() *
            self._rtransform.radius_all()**2 *
            self._int1d.get_weights(self._rtransform.size)
        )

        self._weights = np.ctypeslib.as_ctypes(self.weights)

        self._size = self._rtransform.size

    def integrate(self, f, segments):
        napmo.cext.radial_integrate.restype = c_double
        napmo.cext.radial_integrate.argtypes = [
            POINTER(napmo.RadialGrid),
            c_int,
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')]

        return napmo.cext.radial_integrate(byref(self), segments, f)

    @property
    def size(self):
        return self._size

    @property
    def symbol(self):
        return self._symbol

    @property
    def radii(self):
        return self._radii

    @property
    def rtransform(self):
        return self._rtransform

    @property
    def int1d(self):
        return self._int1d

    def zeros(self):
        return np.zeros(self.shape)
