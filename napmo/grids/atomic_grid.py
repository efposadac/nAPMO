# file: atomic_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.grids.radial_grid import RadialGrid
from napmo.grids.angular_grid import AngularGrid
from napmo.system.c_binding import napmo_library


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

        self._points = np.empty([self.size, 3], dtype=np.float64)
        self._weights = np.empty(self.size, dtype=np.float64)

        offset = 0
        for i in range(n_radial):
            self._points[
                offset:offset + n_angular] = self.radial_grid.points[i] * self.angular_grid.points
            self._weights[
                offset:offset + n_angular] = self.radial_grid.weights[i] * self.angular_grid.weights
            offset += n_angular

        self._points += origin

    def spherical_expansion(self, lmax, f):
        lsize = (lmax + 1) * (lmax + 1)
        output = np.empty([self.radial_grid.size, lsize], dtype=np.float64)
        napmo_library.lebedev_spherical_expansion(
            self.angular_grid.lorder, lmax, self.radial_grid.size, f, output)
        return output

    def evaluate_expansion(self, lmax, expansion):
        output = np.empty(self.size, dtype=np.float64)
        napmo_library.grid_evaluate_atomic_expansion(
            lmax, self.angular_grid.lorder, self.radial_grid.size,
            expansion, self.points, output)
        return output

    def integrate(self, *args):
        args += (np.dstack([self.radial_grid.points**2] *
                           self.angular_grid.lorder).flatten(), )
        f = np.concatenate(args)
        return napmo_library.grid_atomic_integrate(self.size, len(args), self.radial_grid.radii, f, self.weights)

    @property
    def points(self):
        return self._points

    @property
    def weights(self):
        return self._weights

    @property
    def size(self):
        return self._size

    @property
    def origin(self):
        return self._origin


array_1d_double = npct.ndpointer(
    dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(
    dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.lebedev_spherical_expansion.restype = None
napmo_library.lebedev_spherical_expansion.argtypes = [
    c_int,
    c_int,
    c_int,
    array_1d_double,
    array_2d_double
]

napmo_library.grid_evaluate_atomic_expansion.restype = None
napmo_library.grid_evaluate_atomic_expansion.argtypes = [
    c_int,
    c_int,
    c_int,
    array_2d_double,
    array_2d_double,
    array_1d_double
]

napmo_library.grid_atomic_integrate.restype = c_double
napmo_library.grid_atomic_integrate.argtypes = [
    c_int,
    c_int,
    c_double,
    array_1d_double,
    array_1d_double
]
