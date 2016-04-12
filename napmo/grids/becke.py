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

from napmo.system.cext import napmo_library as nl

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.atomic import AtomicGrid
from napmo.utilities.cell import Cell


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
        ("_ncenter", c_int), ("_size", c_int), ("_radii", POINTER(c_double)),
        ("_origin", POINTER(c_double * 3)), ("_points", POINTER(c_double * 3)),
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
            self.atgrids.append(AtomicGrid(n_radial, n_angular, centers[i].get(
                'origin'), centers[i].get('symbol')))

            self.origin[i] = self.atgrids[-1].origin

            self.radii[i] = self.atgrids[-1].radial_grid.radii

            self.points[offset:offset + self.atgrids[-1].size] = self.atgrids[
                -1].points

            self.weights[offset:offset + self.atgrids[-1].size] = self.atgrids[
                -1].weights

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

        nl.becke_weights.restype = None
        nl.becke_weights.argtypes = [
            POINTER(BeckeGrid),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')]

        nl.becke_weights(byref(self), P)
        return P

    def evaluate_decomposition(self, atom, cubic_splines, output, cell=None):
        """
        Evaluate the spherical decomposition for a given atom in the molecular grid:

        :math:`f(\\theta, \\varphi) = \sum_{\ell=0}^{\ell_{max}} \sum_{m=-\ell}^\ell f_{\ell m} \, Y_{\ell m}(\\theta, \\varphi)`

        Args:
            atom (int): Index of the atom.
            cubic_splines (list): List of CubicSpline objects with the spherical expansion.
            output (ndarray): array with the values of the approximated :math:`f(x)`.

        """
        if cell is None:
            cell = Cell(None)

        c_splines = [cubic_splines[i]._this for i in range(len(cubic_splines))]
        splines = np.array(c_splines, dtype=c_void_p)

        nl.eval_decomposition_grid.restype = None
        nl.eval_decomposition_grid.argtypes = [
            npct.ndpointer(c_void_p, flags="C_CONTIGUOUS"),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS'),
            npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS'),
            c_void_p, c_long, c_long]

        nl.eval_decomposition_grid(splines, self.origin[atom, :], output, self.points,
                                   cell._this, len(cubic_splines), self.size)

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
                f[offset:offset + self.atgrids[i].size],
                self.becke_weights[offset:offset + self.atgrids[i].size])
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
