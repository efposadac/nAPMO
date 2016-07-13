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


class PrimitiveGaussian(Structure):

    """
    Defines a Cartesian primitive Gaussian type orbital (GTO). (dict)

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
    _fields_ = [
        ("_l", c_int * 3),
        ("_origin", c_double * 3),
        ("exponent", c_double),
        ("coefficient", c_double),
        ("normalization", c_double)
    ]

    def __init__(self,
                 exponent=0.0,
                 coefficient=1.0,
                 l=np.zeros(3, dtype=np.int32),
                 origin=np.zeros(3, dtype=np.float64)):

        super(PrimitiveGaussian, self).__init__()
        self.exponent = c_double(exponent)
        self.coefficient = c_double(coefficient)
        self._origin[:3] = origin[:]
        self._l[:3] = l[:]
        self.normalization = self.normalize()

    @property
    def origin(self):
        return np.array(self._origin[:3])

    @property
    def l(self):
        return np.array(self._l[:3])

    def normalize(self):
        """
        Calculates the normalization constant of this primitive.
        """
        return napmo.cext.gto_normalize_primitive(byref(self))

    def compute(self, coord):
        """
        Computes the value of the object at ``coord``.
        """
        n_coord = coord.shape[0]
        output = np.empty(n_coord)
        napmo.cext.gto_compute_primitive(
            byref(self), coord, output, n_coord)

        return output

    def overlap(self, other):
        """
        Calculates analytically the overlap integral between primitives.

        Args:
            other (PrimitiveGaussian) : function to perform :math:`<\phi_{self}
             | \phi_{other}>`
        """
        return napmo.cext.gto_overlap_primitive(byref(self), byref(other))

    def _show_compact(self):
        """
        Prints information about the object.
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
