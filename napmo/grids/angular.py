# file: angular_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

import napmo


class AngularGrid(Structure):

    """
    Defines an angular grid based on Lebedev's quadratures.
    """
    _fields_ = [
        ("_lorder", c_int),
        ("_points", POINTER(c_double * 3)),
        ("_weights", POINTER(c_double)),
    ]

    def __init__(self, lorder, spherical=False):
        super(AngularGrid, self).__init__()
        self._lorder = lorder
        self._spherical = spherical

        self.points = np.empty([self.lorder, 3], dtype=np.float64)
        self._points = np.ctypeslib.as_ctypes(self.points)

        self.weights = np.empty(self.lorder, dtype=np.float64)
        self._weights = np.ctypeslib.as_ctypes(self.weights)

        napmo.cext.angular_cartesian.restype = None
        napmo.cext.angular_cartesian.argtypes = [POINTER(AngularGrid)]
        napmo.cext.angular_cartesian(byref(self))

        if spherical:
            napmo.cext.angular_to_spherical.restype = None
            napmo.cext.angular_to_spherical.argtypes = [POINTER(AngularGrid)]
            napmo.cext.angular_to_spherical(byref(self))

    def integrate(self, *args):
        """
        Performs integration over unit sphere by using Lebedev's quadrature.

        Args:
            *args (ndarray): Arrays with the values of function :math:`F(x)` calculated in each point of the sphere.

        Returns:
            integral (double): Integral value.
        """

        napmo.cext.angular_integrate.restype = c_double
        napmo.cext.angular_integrate.argtypes = [
            POINTER(AngularGrid), c_int, npct.ndpointer(
                dtype=np.double, ndim=1, flags='CONTIGUOUS')
        ]

        f = np.concatenate(args)
        return napmo.cext.angular_integrate(byref(self), len(args), f)

    @property
    def lorder(self):
        return self._lorder

    @property
    def spherical(self):
        return self._spherical
