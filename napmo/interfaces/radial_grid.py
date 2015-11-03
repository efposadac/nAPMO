# file: radial_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.utilities.databases import AtomicElementsDatabase
from napmo.interfaces.c_binding import napmo_library


class RadialGrid(object):
    """
    Defines a radial grid for multi-center molecular integration.

    This grid is based on Becke's paper, the transformation requires the covalent radius `rm` of a given atom, such that;

    :math:`r = rm \\frac{1+x}{1-x}`

    References:
        Becke, A. D. A multi-center numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
        size (int): Number of radial points.
        atomic_symbol (str): Atomic symbol for which the grid will be calculated.
    """

    def __init__(self, size, atomic_symbol):
        super(RadialGrid, self).__init__()
        self._rm = AtomicElementsDatabase()[atomic_symbol]['atomic_radii_2']

        self._size = size
        self._symbol = atomic_symbol

        self._points = np.empty(size)
        self._weights = np.empty(size)

        napmo_library.gaussChebyshev(
            size, self._rm, self._points, self._weights)

        self._get_z()
        self._deriv_z()
        self._deriv2_z()

    def _get_z(self):
        """
        Returns the radial points mapped uniformly in the interval [0,1], see Becke's paper.
        """
        self._z = np.empty(self.size)
        napmo_library.gaussChebyshev_get_z(
            self.size, self._rm, self._points, self._z)
        self._step = self._z[0]

    def _deriv_z(self):
        """
        Returns the first derivative of the uniform z grid.
        """
        self._dz = np.empty(self.size)
        napmo_library.gaussChebyshev_deriv_z(
            self.size, self._rm, self._points, self._dz)

    def _deriv2_z(self):
        """
        Returns the second derivative of the uniform z grid.
        """
        self._d2z = np.empty(self.size)
        napmo_library.gaussChebyshev_deriv2_z(
            self.size, self._rm, self._points, self._d2z)

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
    def symbol(self):
        return self._symbol

    @property
    def z(self):
        return self._z

    @property
    def dz(self):
        return self._dz

    @property
    def d2z(self):
        return self._d2z

array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')

napmo_library.gaussChebyshev.restype = None
napmo_library.gaussChebyshev.argtypes = [
    c_int,
    c_double,
    array_1d_double,
    array_1d_double
]

napmo_library.gaussChebyshev_get_z.restype = None
napmo_library.gaussChebyshev_get_z.argtypes = [
    c_int,
    c_double,
    array_1d_double,
    array_1d_double
]

napmo_library.gaussChebyshev_deriv_z.restype = None
napmo_library.gaussChebyshev_deriv_z.argtypes = [
    c_int,
    c_double,
    array_1d_double,
    array_1d_double
]

napmo_library.gaussChebyshev_deriv2_z.restype = None
napmo_library.gaussChebyshev_deriv2_z.argtypes = [
    c_int,
    c_double,
    array_1d_double,
    array_1d_double
]
