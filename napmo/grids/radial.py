# file: radial_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np
from ctypes import *

import napmo


class RadialGrid(object):

    """
    An integration grid for the radial component of a spherical coordinate system
    """

    def __init__(self, size=2, atomic_symbol=None, rtransform=None):

        self._symbol = '--'
        radii = 1.0

        if atomic_symbol in napmo.AtomicElementsDatabase():

            self._symbol = atomic_symbol
            radii = napmo.AtomicElementsDatabase().get(atomic_symbol, {}).get(
                'atomic_radii', 1.0)

        self._rtransform = rtransform

        if rtransform is None:
            self._rtransform = napmo.ChebyshevRadialTransform(radii, size)

        self._this = napmo.cext.RadialGrid_new(self._rtransform._this, self._rtransform.radii)

    def integrate(self, f, segments):
        """
        Perform integration on the radial grid.
        """
        return napmo.cext.RadialGrid_integrate(self._this, segments, f)

    @property
    def size(self):
        return napmo.cext.RadialGrid_get_size(self._this)

    @property
    def symbol(self):
        return self._symbol

    @property
    def radii(self):
        return napmo.cext.RadialGrid_get_radii(self._this)

    @property
    def rtransform(self):
        return self._rtransform

    @property
    def points(self):
        ptr = napmo.cext.RadialGrid_get_points(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

    @property
    def weights(self):
        ptr = napmo.cext.RadialGrid_get_weights(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

    def zeros(self):
        return np.zeros(self.size)
