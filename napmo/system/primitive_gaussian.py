# file: primitive_gaussian.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co
from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


class PrimitiveGaussian(object):

    """
    Defines a Cartesian primitive Gaussian type orbital (GTO).

    Following Obara and Saika (1986) we write an unnormalized primitive
    Cartesian Gaussian function centered at :math:`\\bf A` as

    :math:`\phi ({\\bf r}; \zeta, {\\bf n}, {\\bf A}) = (x - A_x)^{n_x}
    (y - A_y)^{n_y} (z - A_z)^{n_z} \\times \exp[-\zeta({\\bf r}-{\\bf A})^2]`

    where :math:`{\\bf r}` is the coordinate vector of the electron,
    :math:`\zeta` is the orbital exponent, and :math:`\\bf n` is a set of
    non-negative integers. The sum of :math:`n_x`, :math:`n_y`, and :math:`n_z`
    is denoted as :math:`\\bf n` and be referred to as the angular momentum or
    orbital quantum number of the Gaussian function.

    Args:
        exponent (double): GTO exponent.
        coefficient (double): GTO coefficients.
        origin (ndarray) : coordinates (cartesian)
        l (ndarray) : :math:`\\bf n`. Angular moment (x, y, and z components)
    """

    def __init__(self,
                 exponent=0.0,
                 coefficient=1.0,
                 l=np.zeros(3, dtype=np.int32),
                 origin=np.zeros(3, dtype=np.float64)):

        super(PrimitiveGaussian, self).__init__()
        self._this = napmo.cext.PrimitiveGaussian_new(
            l, origin, exponent, coefficient)

    def compute(self, coord):
        """
        Computes the value of the object at ``coord``.
        """
        n_coord = coord.shape[0]
        output = np.empty(n_coord)
        napmo.cext.PrimitiveGaussian_compute(
            self._this, coord, output, n_coord)

        return output

    def overlap(self, other):
        """
        Calculates analytically the overlap integral between two primitives.

        Args:
            other (PrimitiveGaussian) : function to perform :math:`<\phi_{self}
             | \phi_{other}>`
        """
        return napmo.cext.PrimitiveGaussian_overlap(self._this, other._this)

    @property
    def origin(self):
        """
        The center of the function
        """
        output = np.zeros(3)
        napmo.cext.PrimitiveGaussian_get_origin(self._this, output)
        return output

    @property
    def l(self):
        """
        The angular moment of the object
        """
        output = np.zeros(3, dtype=np.int32)
        napmo.cext.PrimitiveGaussian_get_l(self._this, output)
        return output

    @property
    def coefficient(self):
        return napmo.cext.PrimitiveGaussian_get_coeff(self._this)

    @property
    def exponent(self):
        return napmo.cext.PrimitiveGaussian_get_zeta(self._this)

    @property
    def normalization(self):
        return napmo.cext.PrimitiveGaussian_get_norma(self._this)

    def _show_compact(self):
        """
        Prints information about the object in a compact way.
        """
        lvalue = {0: "s", 1: "p", 2: "d", 3: "f", 4: "g"}
        out = '  {0:3s} {1:10.5f} {2:10.5f} {2:10.5f}\n'.format(
            lvalue[sum(self.l)],
            self.exponent,
            self.coefficient,
            self.normalization
        )
        return out

    def __repr__(self):

        out = """
==================================================
Object: {0:9s}
--------------------------------------------------
Origin: {1:<10.5f} {2:<10.5f} {3:<10.5f}
l:      {4:<3d} {5:<3d} {6:<3d}
zeta:   {7:<10.5f}
coeff:  {8:<10.5f}
norma:  {9:<10.5f}
--------------------------------------------------""".format(
            type(self).__name__,
            self.origin[0],
            self.origin[1],
            self.origin[2],
            self.l[0],
            self.l[1],
            self.l[2],
            self.exponent,
            self.coefficient,
            self.normalization,
        )

        return out
