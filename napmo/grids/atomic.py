# file: atomic_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
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
        output = np.zeros([self.radial_grid.size, lsize], dtype=np.float64)

        segments = np.zeros(self.radial_grid.size, dtype=np.int64)
        segments[:] = self.angular_grid.lorder

        f *= self.weights
        f = f.reshape([1, f.size])

        fpp = (f.__array_interface__['data'][
               0] + np.arange(f.shape[0]) * f.strides[0]).astype(np.uintp)

        napmo.cext.dot_multi_moments(f.size, 1, fpp, self.points,
                                     self.origin, lmax, 4, segments, output, lsize)

        output /= self.integrate(segmented=True).reshape(-1, 1)

        counter = 0

        for l in range(0, lmax + 1):
            for m in range(-l, l + 1):
                # proper norm for spherical harmonics
                output[:, counter] *= np.sqrt(4 * np.pi * (2 * l + 1))
                counter += 1

        return output

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
