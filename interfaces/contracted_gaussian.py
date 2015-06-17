# file: contracted_gaussian.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
import numpy as np
import scipy.misc

from interfaces.stack import Stack
from interfaces.primitive_gaussian import PrimitiveGaussian


class ContractedGaussian(dict):
    """
    Defines a linear combination of Cartesian Gaussian type orbitals GTO (dict).

    A contracted Gaussian function is just a linear combination of primitive Gaussians (also termed primitives)
    centered at the same center :math:`{\\bf A}` and with the same momentum indices  :math:`{\\bf n}`
    but with different exponents :math:`\zeta_i`:

    :math:`\phi ({\\bf r}; {\\bf \zeta}, {\\bf C}, {\\bf n}, {\\bf A}) = (x - A_x)^{n_x} (y - A_y)^{n_y} (z - A_z)^{n_z} \\times
    \sum_{i=1}^M C_i \exp[-\zeta_i ({\\bf r}-{\\bf A})^2]`

    Contracted Gaussians form shells the same way as primitives. The contraction coefficients :math:`\\bf C`
    already include normalization constants so that the resulting combination is properly normalized.
    Published contraction coefficients :math:`\\bf c` are linear coefficients for normalized primitives,
    hence the normalization-including contraction coefficients :math:`\\bf C` have to be computed from them as

    :math:`C_i = c_i N(\zeta_i,{\\bf n})`

    where :math:`N` is:

    :math:`N = \\dfrac{1}{\sqrt{<\phi | \phi>}}`

    Args:
        exponents (numpy.ndarray): GTO exponent.
        coefficients (numpy.ndarray): GTO coefficients.
        origin (numpy.ndarray(3)) : coordinates (cartesian)
        l (numpy.ndarray(3)) : :math:`\\bf n`. Angular moment (x, y, and z components)
    """
    def __init__(self, exponents=np.array([0.5]), coefficients=np.array([1.0]), origin=np.array([0.0, 0.0, 0.0]), l=np.array([0, 0, 0])):
        super(ContractedGaussian, self).__init__()

        self["length"] = len(exponents)
        self["l"] = l
        self["origin"] = origin
        self["primitive"] = Stack()

        for (exponent, coefficient) in zip(exponents, coefficients):
            self.get("primitive").push(PrimitiveGaussian(
                exponent, coefficient, l, origin
                ))
        self["normalization"] = 1.0
        aux = self.normalize()
        self["normalization"] = aux

    def normalize(self):
        """
        Normalizes the contraction
        """
        return 1.0 / np.sqrt(self.overlap(self))

    def overlap(self, other):
        """
        Calculates the overlap integral between two contractions.

        Args:
            other (ContractedGaussian) : Contracted function to perform :math:`<\phi_{self} | \phi_{other}>`
        """
        output = 0.0
        for pa in self.get('primitive'):
            for pb in other.get('primitive'):
                output += pa.overlap(pb)

        output *= self.get('normalization') * other.get('normalization')

        return output

    def compute(self, coord):
        """
        Computes the value of the contracted Gaussian at ``coord``.
        """
        RP = coord - self.get('origin')
        RP2 = RP.dot(RP)

        factor = 1.0
        for i in range(3):
            factor *= RP[i]**self.get('l')[i]

        output = 0.0
        for primitive in self.get('primitive'):
            output += (
                primitive.get('coefficient') *
                primitive.get('normalization') *
                self.get('normalization') *
                factor *
                np.exp(-primitive.get('exponent') * RP2)
            )
        return output

    def show(self):
        """
        Prints the contents of the object.
        """
        print("  origin: ", self.get('origin'))
        print("  length: ", self.get('length'))
        print("  normalization: ", self.get('normalization'))
        print("  *** Primitives info: ")
        i = 1
        for primitive in self.get('primitive'):
            print(i)
            primitive.show()
            i += 1
