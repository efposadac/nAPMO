# file: atomic_grid.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

import napmo


class AtomicGrid(object):

    """
    Defines a spherical grid for each center in the system.
    """

    def __init__(self, nrad, nang, origin, atomic_symbol, rtransform=None):
        super(AtomicGrid, self).__init__()

        self.radial_grid = napmo.RadialGrid(
            nrad, atomic_symbol, rtransform=rtransform)

        self.angular_grid = napmo.AngularGrid(nang)

        self._symbol = atomic_symbol

        self._this = napmo.cext.AtomicGrid_new(
            self.angular_grid._this, self.radial_grid._this, origin)

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

        if lmax < 0:
            raise ValueError('lmax can not be negative.')

        lsize = (lmax + 1) * (lmax + 1)

        ptr = napmo.cext.AtomicGrid_spherical_expansion(self._this, f, lmax)

        return np.ctypeslib.as_array(ptr, shape=(self.radial_grid.size, lsize))

    def evaluate_expansion(self, lmax, expansion):
        """
        Evaluate the spherical expansion for one atom (do not use with molecules):

        :math:`f(\\theta, \\varphi) = \sum_{\ell=0}^{\ell_{max}} \sum_{m=-\ell}^\ell f_{\ell m} \, Y_{\ell m}(\\theta, \\varphi)`

        Args:
            lmax (int): Maximum :math:`\ell` order of the expansion.
            expansion (ndarray): Spherical expansion array with shape (nrad, lsize), where :math:`\ell_{size} = (\ell_{max} + 1)^2`

        Returns:
            output (ndarray): array with the values of the approximated :math:`f(x)`.
        """
        output = np.empty(self.size, dtype=np.float64)

        napmo.cext.AngularGrid_eval_expansion(
            self.angular_grid._this, lmax, self.radial_grid.size, expansion, output)

        return output

    def spherical_average(self, *args):
        """
        Computes the spherical average of the product of given functions

        Returns:
            output (array) : the spherical average
        """

        f = np.concatenate(args)
        s = self.integrate(segmented=True)
        output = self.integrate(f, segmented=True)
        output /= s

        return output

    def integrate(self, *args, **kwargs):
        """
        Perform an integration of function :math:`f` using AtomicGrid.

        Args:
            f (ndarray): array of :math:`f` computed in all grid points.
            segmented (logical): whether to calculate the integral in a segmented way along the number of radial points.

        Returns:
            integral: float or array in the segmented case.
        """
        segmented = kwargs.pop('segmented', False)

        args += (self.weights,)
        nfunc = len(args)

        if segmented:
            nseg = self.radial_grid.size
            sseg = self.angular_grid.lorder
        else:
            nseg = 1
            sseg = self.size

        f = np.concatenate(args)

        ptr = napmo.cext.AtomicGrid_integrate(
            self._this, nfunc, nseg, sseg, f)

        return np.ctypeslib.as_array(ptr, shape=(nseg,))

    @property
    def points(self):
        ptr = napmo.cext.AtomicGrid_get_points(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size, 3))

    @property
    def weights(self):
        ptr = napmo.cext.AtomicGrid_get_weights(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

    @property
    def origin(self):
        ptr = napmo.cext.AtomicGrid_get_origin(self._this)
        return np.ctypeslib.as_array(ptr, shape=(3,))

    @property
    def size(self):
        return napmo.cext.AtomicGrid_get_size(self._this)

    @property
    def radii(self):
        return napmo.cext.AtomicGrid_get_radii(self._this)

    @property
    def symbol(self):
        """
        Owner of the grid
        """
        return self._symbol
