# file: becke.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *
from scipy.io import FortranFile

import napmo


class BeckeGrid(object):

    """
    This class creates the Becke grid.

    References:
        Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
        species (dict): Species to be calculated. see MolecularSystem
        n_radial (int, optional): Number of radial points. Default is 40
        n_angular (int, optional): Number of angular points. Default is 110
    """

    def __init__(self, species, n_radial=40, n_angular=110, rtransform=None, file=None, ablmax=None, abldep=None):
        super(BeckeGrid, self).__init__()

        assert isinstance(species, dict)

        centers = species.get('particles')
        ncenter = len(centers)
        self._symbol = species.get('symbol')

        if file is not None:
            f = FortranFile(file, 'r')
            size = f.read_ints(dtype=np.int64)[0]
            pw = f.read_reals(dtype=np.float64)
            pw = pw.reshape(int(size / 4), 4, order='F')
            pw = np.asanyarray(pw, order='C')
            size = pw.shape[0]

            self._this = napmo.cext.BeckeGrid_from_points(pw, size, ncenter)
        else:

            if ablmax is None:
                ablmax = 0

            if abldep is None:
                abldep = 1.0e-6

            self._nrad = n_radial
            self._nang = n_angular

            self.atgrids = [napmo.AtomicGrid(n_radial, n_angular, center.get(
                'origin'), center.get('symbol'), rtransform=rtransform)
                for center in centers]

            atgrids_ptr = np.array(
                [atgrid._this for atgrid in self.atgrids], dtype=c_void_p)

            self._this = napmo.cext.BeckeGrid_new(atgrids_ptr, ncenter, ablmax, abldep)

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
            cell = napmo.Cell(None)

        c_splines = [cs._this for cs in cubic_splines]
        splines = np.array(c_splines, dtype=c_void_p)

        napmo.cext.eval_decomposition_grid(splines, self.origin[atom, :], output, self.points,
                                           cell._this, len(cubic_splines), self.size)

    def integrate(self, f):
        """
        Perform an integration of function :math:`F(r)` using BeckeGrid.

        Args:
            f (ndarray): array of F computed in all grid points.
        """
        return napmo.cext.BeckeGrid_integrate(self._this, f)

    def show(self):
        """
        Prints information of the object.
        """
        print("\nGrid Information:", self._symbol)
        print("-" * (18 + len(self._symbol)))
        print("Centers: ", self.ncenter)
        print("Size: ", self.size)

    @property
    def becke_weights(self):
        """Computes the Becke weights :math:`w(r)` for the entire grid as described in eq. 22 Becke, 1988.

        References:
            Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

        Returns:
            P (ndarray): The value of cell_function (eq. 13, Becke, 1988)
        """
        ptr = napmo.cext.BeckeGrid_get_becke_weights(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

    @property
    def points(self):
        ptr = napmo.cext.BeckeGrid_get_points(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size, 3))

    @property
    def weights(self):
        ptr = napmo.cext.BeckeGrid_get_weights(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

    @property
    def origin(self):
        ptr = napmo.cext.BeckeGrid_get_origin(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.ncenter, 3))

    @property
    def size(self):
        return napmo.cext.BeckeGrid_get_size(self._this)

    @property
    def radii(self):
        ptr = napmo.cext.BeckeGrid_get_radii(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.ncenter,))

    @property
    def ncenter(self):
        """
        Number of particles of a given species in the molecular grid
        """
        return napmo.cext.BeckeGrid_get_ncenter(self._this)

    @property
    def nrad(self):
        return self._nrad

    @property
    def nang(self):
        return self._nang

    @property
    def symbol(self):
        return self._symbol
