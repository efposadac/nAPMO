# file: primitive_gaussian.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
from __future__ import print_function

import numpy as np
import numpy.ctypeslib as npct
from ctypes import *

from napmo.system.cext import napmo_library


class PrimitiveGaussian(Structure):
    """
    Defines a Cartesian primitive Gaussian type orbital (GTO). (dict)

    Following Obara and Saika (1986) we write an unnormalized primitive Cartesian Gaussian function centered at :math:`\\bf A` as

    :math:`\phi ({\\bf r}; \zeta, {\\bf n}, {\\bf A}) = (x - A_x)^{n_x} (y - A_y)^{n_y} (z - A_z)^{n_z}
    \\times \exp[-\zeta({\\bf r}-{\\bf A})^2]`

    where :math:`{\\bf r}` is the coordinate vector of the electron, :math:`\zeta` is the orbital exponent, and :math:`\\bf n`
    is a set of non-negative integers. The sum of :math:`n_x`, :math:`n_y`, and :math:`n_z` is denoted as :math:`\\bf n`
    and be referred to as the angular momentum or orbital quantum number of the Gaussian function.

    Args:
        exponent (float64): GTO exponent.
        coefficient (float64): GTO coefficients.
        origin (numpy.ndarray(3)) : coordinates (cartesian)
        l (numpy.ndarray(3)) : :math:`\\bf n`. Angular moment (x, y, and z components)
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
        return napmo_library.gto_normalize_primitive(byref(self))

    def compute(self, coord):
        """
        Computes the value of the object at ``coord``.
        """
        n_coord = coord.shape[0]
        output = np.empty(n_coord)
        napmo_library.gto_compute_primitive(
            byref(self), coord, output, n_coord)

        return output

    def overlap(self, other):
        """
        Calculates analytically the overlap integral between primitives.

        Args:
            other (PrimitiveGaussian) : function to perform :math:`<\phi_{self} | \phi_{other}>`
        """
        return napmo_library.gto_overlap_primitive(byref(self), byref(other))

    def show(self):
        """
        Prints information about the object.
        """
        print("    origin: ", self.origin)
        print("    exponent: ", self.exponent)
        print("    coefficient: ", self.coefficient)
        print("    angular moment: ", self.l)
        print("    normalization: ", self.normalization)


array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.gto_normalize_primitive.restype = c_double
napmo_library.gto_normalize_primitive.argtypes = [
    POINTER(PrimitiveGaussian)
]

napmo_library.gto_compute_primitive.restype = None
napmo_library.gto_compute_primitive.argtypes = [
    POINTER(PrimitiveGaussian),
    array_2d_double, array_1d_double,
    c_int
]

napmo_library.gto_overlap_primitive.restype = c_double
napmo_library.gto_overlap_primitive.argtypes = [
    POINTER(PrimitiveGaussian),
    POINTER(PrimitiveGaussian)
]
