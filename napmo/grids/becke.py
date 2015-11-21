# file: becke.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.system.cext import napmo_library
from napmo.system.molecular_system import MolecularSystem
from napmo.grids.atomic import AtomicGrid


class BeckeGrid(Structure):
    """
    This class creates the Becke grid.

    References:
        Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
        molecule (MolecularSystem): Molecular system to be calculated.
        n_radial (int, optional): Number of radial points. Default is 40
        n_angular (int, optional): Number of angular points. Default is 110
    """
    _fields_ = [
        ("_ncenter", c_int),
        ("_size", c_int),
        ("_radii", POINTER(c_double)),
        ("_origin", POINTER(c_double * 3)),
        ("_points", POINTER(c_double * 3)),
        ("_weights", POINTER(c_double))
    ]

    def __init__(self, molecule, n_radial=40, n_angular=110):
        super(BeckeGrid, self).__init__()

        assert isinstance(molecule, MolecularSystem)

        centers = molecule.get('atoms')

        self._ncenter = len(centers)
        self._size = n_radial * n_angular * self.ncenter

        self.radii = np.empty(self.ncenter, dtype=np.float64)
        self._radii = np.ctypeslib.as_ctypes(self.radii)

        self.origin = np.empty([self.ncenter, 3], dtype=np.float64)
        self._origin = np.ctypeslib.as_ctypes(self.origin)

        self.points = np.empty([self.size, 3], dtype=np.float64)
        self._points = np.ctypeslib.as_ctypes(self.points)

        self.weights = np.empty(self.size, dtype=np.float64)
        self._weights = np.ctypeslib.as_ctypes(self.weights)

        self.atgrids = []

        offset = 0
        for i in range(self.ncenter):
            self.atgrids.append(
                AtomicGrid(n_radial, n_angular, centers[i].get('origin'),
                           centers[i].get('symbol')))

            self.origin[i] = self.atgrids[-1].origin

            self.radii[i] = self.atgrids[-1].radial_grid.radii

            self.points[offset:offset +
                        self.atgrids[-1].size] = self.atgrids[-1].points

            self.weights[offset:offset +
                         self.atgrids[-1].size] = self.atgrids[-1].weights

            offset += self.atgrids[-1].size

        self.becke_weights = self._becke_weights()

    def _becke_weights(self):
        """Computes the Becke weights :math:`w(r)` for the entire grid as described in eq. 22 Becke, 1988.

        References:
            Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

        Returns:
            P (numpy.darray): The value of cell_function (eq. 13, Becke, 1988)
        """
        P = np.empty(self.size, dtype=np.float64)
        napmo_library.becke_weights(byref(self), P)
        return P

    def integrate(self, f):
        """
        Perform an integration of function :math:`F(r)` using BeckeGrid.

        Args:
            f (ndarray): array of F computed in all grid points.
        """
        offset = 0
        integral = 0.0
        for i in range(self.ncenter):
            integral += self.atgrids[i].integrate(
                f[offset:offset + self.atgrids[i].size], self.becke_weights[offset:offset + self.atgrids[i].size])
            offset += self.atgrids[i].size
        return integral

    def show(self):
        """
        Prints information of the object.
        """
        print("Grid Information:")
        print("-----------------")
        print("Centers: ", self.ncenter)
        print("Size: ", self.size)

    @property
    def ncenter(self):
        return self._ncenter

    @property
    def size(self):
        return self._size


array_1d_double = npct.ndpointer(
    dtype=np.double, ndim=1, flags='CONTIGUOUS')

napmo_library.becke_weights.restype = None
napmo_library.becke_weights.argtypes = [
    POINTER(BeckeGrid),
    array_1d_double
]
