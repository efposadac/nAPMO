# file: radial_grid.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np
from ctypes import *

import napmo


class RadialGrid(object):

    """
    An integration grid for the radial component of a spherical coordinate system
    """

    def __init__(self, size=2, atomic_symbol="--", rtransform=None):

        self._symbol = atomic_symbol

        # Hydrogen by default
        radii = napmo.AtomicElementsDatabase().get("H", {}).get(
            'atomic_radii_2', 1.0)

        # Isotopes case
        aux = atomic_symbol.find("_")
        if aux > 0:
            atom = atomic_symbol[0:aux]
            radii = napmo.AtomicElementsDatabase().get(atom, {}).get(
                'atomic_radii', 1.0)

        # Atoms case
        if atomic_symbol in napmo.AtomicElementsDatabase():

            radii = napmo.AtomicElementsDatabase().get(atomic_symbol, {}).get(
                'atomic_radii', 1.0)

        self._rtransform = rtransform

        print("Radii used for " + atomic_symbol + ": ", radii)

        if rtransform is None:
            self._rtransform = napmo.ChebyshevRadialTransform(radii, size)

        self._this = napmo.cext.RadialGrid_new(self._rtransform._this, self._rtransform.radii)

    def integrate(self, f, segments):
        """
        Perform integration on the radial grid.
        """
        return napmo.cext.RadialGrid_integrate(self._this, segments, f)

    def deriv2(self, f):
        """
        Perform integration on the radial grid.
        """
        ptr = napmo.cext.RadialGrid_deriv2(self._this, f)
        return np.ctypeslib.as_array(ptr, shape=(self.size,))

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
