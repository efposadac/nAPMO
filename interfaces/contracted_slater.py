# file: contracted_slater.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it
from __future__ import division
import numpy as np

from interfaces.stack import Stack
from interfaces.primitive_slater import PrimitiveSlater


class ContractedSlater(dict):
    """
    Defines a Slater Type Orbital (STO)  as a combination of PrimitiveSlater Functions (spherical coordinates)
    """
    def __init__(self, exponents=np.array([0.5]), coefficients=np.array([1.0]), origin=np.array([0.0, 0.0, 0.0]), n=np.array([0]), l=0, m=0):
        super(ContractedSlater, self).__init__()

        self["length"] = len(exponents)
        self["origin"] = origin
        self["primitive"] = Stack()
        self['n'] = n
        self['l'] = l
        self['m'] = m

        for (exponent, coefficient, ni) in zip(exponents, coefficients, n):
            self.get("primitive").push(PrimitiveSlater(
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
        """
        RP = np.zeros(3)
        RP = coord - self.get('origin')

        xy = RP[0]**2 + RP[1]**2
        r = np.sqrt(xy + RP[2]**2)
        theta = np.arctan2(np.sqrt(xy), RP[2])
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
