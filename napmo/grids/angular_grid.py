# file: angular_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.utilities.angular_quadratures import lebedev
from napmo.system.c_binding import napmo_library


class AngularGrid(object):
    """
    Defines an angular grid based on Lebedev's quadratures.
    """

    def __init__(self, lorder, spherical=False):
        super(AngularGrid, self).__init__()
        self._lorder = lorder
        self._spherical = spherical
        self._points, self._weights = lebedev(lorder, spherical)

    def integrate(self, *args):
        f = np.concatenate(args)
        return napmo_library.lebedev_integrate(self.lorder, len(args), f, self.weights)

    @property
    def points(self):
        return self._points

    @property
    def weights(self):
        return self._weights

    @property
    def lorder(self):
        return self._lorder

    @property
    def spherical(self):
        return self._spherical

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')

napmo_library.lebedev_integrate.restype = c_double
napmo_library.lebedev_integrate.argtypes = [
    c_int,
    c_int,
    array_1d_double,
    array_1d_double
]
