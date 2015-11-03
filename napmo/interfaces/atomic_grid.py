# file: atomic_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy.ctypeslib as npct

from napmo.interfaces.radial_grid import *
from napmo.interfaces.angular_grid import *
from napmo.interfaces.c_binding import napmo_library


class AtomicGrid(object):
    """
    docstring for AtomicGrid
    """

    def __init__(self, n_radial, n_angular, origin, atomic_symbol):
        super(AtomicGrid, self).__init__()
        self._size = n_radial * n_angular
        self._symbol = atomic_symbol
        self._origin = origin

        self.radial_grid = RadialGrid(n_radial, atomic_symbol)
        self.angular_grid = AngularGrid(n_angular)

        self._points = np.empty([self.size, 3])
        self._weights = np.empty(self.size)

        offset = 0
        for i in range(n_radial):
            self._points[
                offset:offset + n_angular] = self.radial_grid.points[i] * self.angular_grid.points
            self._weights[
                offset:offset + n_angular] = self.radial_grid.weights[i] * self.angular_grid.weights
            offset += n_angular

    def spherical_expansion(self, lmax, f):
        lsize = (lmax + 1) * (lmax + 1)
        output = np.empty([self.radial_grid.size, lsize])
        napmo_library.lebedev_spherical_expansion(
            self.angular_grid.lorder, lmax, self.radial_grid.size, f, output)
        return output

    def evaluate_expansion(self):
        pass

    def integrate(self):
        pass

    @property
    def points(self):
        return self._points

    @property
    def weights(self):
        return self._weights

    @property
    def size(self):
        return self._size

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.lebedev_spherical_expansion.restype = None
napmo_library.lebedev_spherical_expansion.argtypes = [
    c_int,
    c_int,
    c_int,
    array_1d_double,
    array_2d_double
]
