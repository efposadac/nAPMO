# file: atomic_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.grids.radial import RadialGrid
from napmo.grids.angular import AngularGrid
from napmo.system.cext import napmo_library


class AtomicGrid(Structure):
    """
    Defines a spherical grid for each center in the system.
    """
    _fields_ = [
        ("_size", c_int),
        ("_radii", c_double),
        ("_origin", POINTER(c_double * 3)),
        ("_points", POINTER(c_double * 3)),
        ("_weights", POINTER(c_double)),
    ]

    def __init__(self, nrad, nang, origin, atomic_symbol):
        super(AtomicGrid, self).__init__()
        self._size = nrad * nang

        self.origin = np.array([origin], dtype=np.float64)
        self._origin = np.ctypeslib.as_ctypes(self.origin)

        self.radial_grid = RadialGrid(nrad, atomic_symbol)
        self.angular_grid = AngularGrid(nang)

        self._radii = self.radial_grid.radii

        self.points = np.empty([self.size, 3], dtype=np.float64)
        self._points = np.ctypeslib.as_ctypes(self.points)

        self.weights = np.empty(self.size, dtype=np.float64)
        self._weights = np.ctypeslib.as_ctypes(self.weights)

        self._symbol = atomic_symbol

        offset = 0
        for i in range(nrad):
            self.points[offset:offset + nang] = (
                self.radial_grid.points[i] * self.angular_grid.points
            )
            self.weights[offset:offset + nang] = (
                self.radial_grid.weights[i] * self.angular_grid.weights *
                self.radial_grid.points[i] * self.radial_grid.points[i]
            )
            offset += nang

        self.points += origin

    def spherical_expansion(self, lmax, f):
        """
        Performs the spherical expansion:

        :math:`f_{\ell m}=\int_{\Omega} f(\\theta,\\varphi)\, Y_{\ell m}(\\theta,\\varphi)\,d\Omega`

        Args:
            lmax (int): Maximum :math:`\ell` order of the expansion.
            f (ndarray): Array with the values of function :math:`f(x)` calculated in each point of the atomic grid.

        Returns:
            integral (ndarray): Spherical expansion array with shape (nrad, lsize), where :math:`\ell_{size} = (\ell_{max} + 1)^2`
        """
        lsize = (lmax + 1) * (lmax + 1)
        output = np.empty([self.radial_grid.size, lsize], dtype=np.float64)
        napmo_library.angular_spherical_expansion(
            byref(self.angular_grid), lmax, self.radial_grid.size, f, output)
        return output

    def evaluate_expansion(self, lmax, expansion):
        """
        Evaluate the spherical expansion:

        :math:`f(\\theta, \\varphi) = \sum_{\ell=0}^{\ell_{max}} \sum_{m=-\ell}^\ell f_{\ell m} \, Y_{\ell m}(\\theta, \\varphi)`

        Args:
            lmax (int): Maximum :math:`\ell` order of the expansion.
            expansion (ndarray): Spherical expansion array with shape (nrad, lsize), where :math:`\ell_{size} = (\ell_{max} + 1)^2`

        Returns:
            output (ndarray): array with the values of the approximated :math:`f(x)`.
        """
        output = np.empty(self.size, dtype=np.float64)
        napmo_library.angular_eval_expansion(
            byref(self.angular_grid), lmax, self.radial_grid.size, expansion, output)
        return output

    def integrate(self, *args):
        """
        Perform an integration of function :math:`f` using AtomicGrid.

        Args:
            f (ndarray): array of :math:`f` computed in all grid points.
        """
        f = np.concatenate(args)
        integral = napmo_library.atomic_grid_integrate(
            byref(self), len(args), f)
        return integral

    @property
    def size(self):
        return self._size

    @property
    def symbol(self):
        return self._symbol


array_1d_double = npct.ndpointer(
    dtype=np.double, ndim=1, flags='CONTIGUOUS')

array_2d_double = npct.ndpointer(
    dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.angular_spherical_expansion.restype = None
napmo_library.angular_spherical_expansion.argtypes = [
    POINTER(AngularGrid),
    c_int,
    c_int,
    array_1d_double,
    array_2d_double
]

napmo_library.angular_eval_expansion.restype = None
napmo_library.angular_eval_expansion.argtypes = [
    POINTER(AngularGrid),
    c_int,
    c_int,
    array_2d_double,
    array_1d_double
]

napmo_library.atomic_grid_integrate.restype = c_double
napmo_library.atomic_grid_integrate.argtypes = [
    POINTER(AtomicGrid),
    c_int,
    array_1d_double
]
