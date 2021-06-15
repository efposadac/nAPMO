# file: angular_grid.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

import napmo


class AngularGrid(object):

    """
    Defines an angular grid based on Lebedev's quadratures.
    """

    def __init__(self, lorder, spherical=False):
        super(AngularGrid, self).__init__()

        self._this = napmo.cext.AngularGrid_new(lorder)

        self._spherical = spherical

        if spherical:
            napmo.cext.AngularGrid_spherical(self._this)

    def integrate(self, *args):
        """
        Performs integration over unit sphere by using Lebedev's quadrature.

        Args:
            *args (ndarray): Arrays with the values of function :math:`F(x)` calculated in each point of the sphere.

        Returns:
            integral (double): Integral value.
        """
        f = np.concatenate(args)
        return napmo.cext.AngularGrid_integrate(self._this, len(args), f)

    @property
    def lorder(self):
        """
        Order of the Lebedev's quadrature
        """
        return napmo.cext.AngularGrid_get_lorder(self._this)

    @property
    def points(self):
        ptr = napmo.cext.AngularGrid_get_points(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.lorder, 3))

    @property
    def weights(self):
        ptr = napmo.cext.AngularGrid_get_weights(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.lorder,))

    @property
    def spherical(self):
        """
        Whether the points are in spherical coordinates or not
        """
        return self._spherical
