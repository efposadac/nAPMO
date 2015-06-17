# file: primitive_gaussian.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
import numpy as np
import scipy.misc

from utilities import analytical_integration as aint


class PrimitiveGaussian(dict):
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
    def __init__(self, exponent=0.0, coefficient=1.0, l=np.array([0, 0, 0]), origin=np.array([0.0, 0.0, 0.0])):
        super(PrimitiveGaussian, self).__init__()
        self["exponent"] = exponent
        self["coefficient"] = coefficient
        self["l"] = l
        self["origin"] = np.array(origin)
        self["normalization"] = self.normalize()

    def normalize(self):
        """
        Calculates the normalization constant of this primitive.
        """
        output = ((2.0 * self.get('exponent')/np.pi)**0.75) / np.sqrt(
                    scipy.misc.factorial2(np.abs(2 * self.get('l')[0] - 1)) *
                    scipy.misc.factorial2(np.abs(2 * self.get('l')[1] - 1)) *
                    scipy.misc.factorial2(np.abs(2 * self.get('l')[2] - 1)) /
                    ((4.0 * self.get('exponent'))**np.sum(self.get('l'))))

        return output

    def compute(self, coord):
        """
        Computes the value of the object at ``coord``.
        """
        RP = coord - self.get('origin')
        RP2 = RP.dot(RP)

        factor = 1.0
        for i in range(3):
            factor *= RP[i]**self.get('l')[i]
        output = (
                    self.get('coefficient') *
                    self.get('normalization') *
                    factor *
                    np.exp(-self.get('exponent') * RP2)
                )

        return output

    def overlap(self, other):
        """
        Calculates analytically the overlap integral between primitives.

        Args:
            other (PrimitiveGaussian) : function to perform :math:`<\phi_{self} | \phi_{other}>`
        """
        gamma = self.get('exponent') + other.get('exponent')
        gammaInv = 1.0/gamma

        AB = self.get('origin') - other.get('origin')
        AB2 = AB.dot(AB)

        P0 = np.zeros(3)
        PA = np.zeros(3)
        PB = np.zeros(3)

        P0 = (self.get('exponent') * self.get('origin') + other.get('exponent') * other.get('origin')) * gammaInv
        PA = P0 - self.get('origin')
        PB = P0 - other.get('origin')

        preFactor = np.exp(- self.get('exponent')*other.get('exponent')*AB2*gammaInv) * (
            np.sqrt(np.pi*gammaInv) * np.pi * gammaInv *
            self.get('coefficient') * other.get('coefficient') *
            self.get('normalization') * other.get('normalization')
            )

        # recursion
        x, y, z = aint.obaraSaika_recursion(PA, PB, gamma, np.sum(self.get('l'))+2, np.sum(other.get('l'))+2)

        x0 = x[self.get('l')[0], other.get('l')[0]]
        y0 = y[self.get('l')[1], other.get('l')[1]]
        z0 = z[self.get('l')[2], other.get('l')[2]]

        # Calculating integrals for primitives
        return preFactor*x0*y0*z0

    def show(self):
        """
        Prints information about the object.
        """
        print("    origin: ", self.get('origin'))
        print("    exponent: ", self.get('exponent'))
        print("    coefficient: ", self.get('coefficient'))
        print("    angular moment: ", self.get('l'))
        print("    normalization: ", self.get('normalization'))
