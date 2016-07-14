# file: radial_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np
from ctypes import *
import napmo


class RadialGridCheb(Structure):

    """
    Defines a radial grid for multi-center molecular integration.

    This grid is based on Becke's paper, the transformation requires the covalent radius ``rm`` of a given atom, such that;

    :math:`r = rm \\dfrac{1+x}{1-x}`

    When the object is instanced calculates the mapped uniform grid ``z`` and its first and second derivatives needed to calculate
    coulomb potentials.

    References:
        Becke, A. D. A multi-center numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
        size (int): Number of radial points.
        atomic_symbol (str): Owner of the grid
    """
    _fields_ = [
        ("_size", c_int),
        ("_radii", c_double),
        ("_points", POINTER(c_double)),
        ("_weights", POINTER(c_double)),
        ("_z", POINTER(c_double)),
        ("_dz", POINTER(c_double)),
        ("_d2z", POINTER(c_double))
    ]

    def __init__(self, size, atomic_symbol):
        super(RadialGridCheb, self).__init__()
        self._radii = napmo.AtomicElementsDatabase()[atomic_symbol][
            'atomic_radii_2']
        self._size = size
        self._symbol = atomic_symbol

        self.points = np.empty(self.size, dtype=np.float64)
        self._points = np.ctypeslib.as_ctypes(self.points)

        self.weights = np.empty(self.size, dtype=np.float64)
        self._weights = np.ctypeslib.as_ctypes(self.weights)

        napmo.cext.radial_init(byref(self))

        self._get_z()
        self._deriv_z()
        self._deriv2_z()

    def _get_z(self):
        """
        Returns the radial points mapped uniformly in the interval :mat:`[0,1]`, see Becke's paper.
        """
        self.z = np.empty(self.size, dtype=np.float64)
        self._z = np.ctypeslib.as_ctypes(self.z)

        napmo.cext.radial_get_z(byref(self))

    def _deriv_z(self):
        """
        Returns the first derivative of the uniform ``z`` grid.
        """
        self.dz = np.empty(self.size, dtype=np.float64)
        self._dz = np.ctypeslib.as_ctypes(self.dz)

        napmo.cext.radial_deriv_z(byref(self))

    def _deriv2_z(self):
        """
        Returns the second derivative of the uniform ``z`` grid.
        """
        self.d2z = np.empty(self.size, dtype=np.float64)
        self._d2z = np.ctypeslib.as_ctypes(self.d2z)

        napmo.cext.radial_deriv2_z(byref(self))

    @property
    def size(self):
        """
        The number of points in the grid
        """
        return self._size

    @property
    def symbol(self):
        """
        Particle owner of the grid
        """
        return self._symbol

    @property
    def radii(self):
        """
        The grid points
        """
        return self._radii
