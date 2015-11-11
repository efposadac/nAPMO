# file: contracted_slater.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
from __future__ import print_function

import numpy as np

from napmo.system.primitive_slater import PrimitiveSlater


class ContractedSlater(dict):
    """
    Defines a Slater Type Orbital (STO)  as a combination of PrimitiveSlater Functions (spherical coordinates). (dict)

    Args:
        exponents (ndarray): Slater exponent.
        coefficient (ndarray): Slater coefficients.
        origin (ndarray) : coordinates (spherical)
        n (ndarray): Quantum number n = 1, 2, 3, etc. (one for each primitive)
        l (int): Quantum number l = S=0, P=1, D=2, etc.
        m (int): Quantum number m = order of the harmonic, m <= l
    """

    def __init__(self,
                 exponents=np.array([0.5], dtype=np.float64),
                 coefficients=np.array([1.0], dtype=np.float64),
                 origin=np.zeros(3, dtype=np.float64),
                 n=np.zeros(1, dtype=np.int32),
                 l=0,
                 m=0):

        super(ContractedSlater, self).__init__()

        self["length"] = len(exponents)
        self["origin"] = origin
        self["primitive"] = []
        self['n'] = n
        self['l'] = l
        self['m'] = m

        for (exponent, coefficient, ni) in zip(exponents, coefficients, n):
            self.get("primitive").append(PrimitiveSlater(
                exponent, coefficient, ni, l, m, origin
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
            other (ContractedSlater) : Contracted function to perform :math:`<\phi_{self} | \phi_{other}>`
        """
        output = 0.0
        for pa in self.get('primitive'):
            for pb in other.get('primitive'):
                output += pa.overlap(pb)

        output *= self.get('normalization') * other.get('normalization')

        return output

    def compute(self, coord):
        """
        Computes the value of the contracted STO at ``coord``.

        Args:
            coord (ndarray) : array with the points where the function will be calculated.
        """
        RP = np.zeros(3, dtype=np.float64)
        RP = coord - self.get('origin')

        if coord.ndim > 1:
            xy = RP[:, 0]**2 + RP[:, 1]**2
            r = np.sqrt(xy + RP[:, 2]**2)
            theta = np.arctan2(np.sqrt(xy), RP[:, 2])
            phi = np.arctan2(RP[:, 1], RP[:, 0])

        else:
            xy = RP[0]**2 + RP[1]**2
            r = np.sqrt(xy + RP[2]**2)
            theta = np.arctan2(np.sqrt(xy), RP[2])
            # if RP[0] == 0.0 and RP[1] == 0.0:
            #     phi = np.arctan2(RP[2], RP[1])
            # else:
            #     phi = np.arctan2(RP[1], RP[0])
            phi = np.arctan2(RP[1], RP[0])

        output = 0.0
        for primitive in self.get('primitive'):
            output += (
                primitive.get('normalization') *
                primitive.get('coefficient') *
                primitive.radial(r) *
                primitive.spherical(theta, phi) *
                self.get('normalization')
            )

        return output

    def show(self):
        """
        Prints the contents of the object.
        """
        print("  origin: ", self.get('origin'))
        print("  length: ", self.get('length'))
        print("  n: ", self.get('n'))
        print("  l: ", self.get('l'))
        print("  m: ", self.get('m'))
        print("  normalization: ", self.get('normalization'))
        print("")
        print("  *** Primitives info: ")
        print("")
        i = 1
        for primitive in self.get('primitive'):
            print(i)
            primitive.show()
            i += 1
