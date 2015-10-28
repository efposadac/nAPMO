# file: primitive_slater.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
from __future__ import division
from __future__ import print_function

from scipy.special import sph_harm
import scipy.misc
import numpy as np


class PrimitiveSlater(dict):
    """Slater-Type orbitals. (dict):

    The slater orbital function is:

    :math:`\chi_{p \lambda \\alpha}(r, \\theta, \phi) = R_{ \lambda p}(r) Y_{ \lambda \\alpha}( \\theta, \phi)`

    Where,

    :math:`R_{\lambda p} =
    [(2n_{\lambda p})!]^{-1/2}(2\zeta_{\lambda p})^{n_{\lambda p}+1/2} r^{\lambda p -1} e^{-\zeta_{\lambda p} r}`

    And,

    :math:`Y_{\lambda \\alpha}(\\theta,\phi)` is an spherical harmonic.

    References:
        1. E. Clementi, C. Roetti, Roothaan-Hartree-Fock atomic wavefunctions. \
At. Data Nucl. Data Tables. 14. (3)-(4), 177-478 (1974).

    Args:
        zeta (float64): Slater exponent.
        n (int): Quantum number n = 1, 2, 3, etc.
        l (int): Quantum number l = S=0, P=1, D=2, etc.
        m (int): Quantum number m = order of the harmonic, m <= l
        origin (numpy.ndarray(3)) : coordinates (spherical)

    """
    def __init__(self, exponent, coefficient=1.0, n=1, l=0, m=0, origin=np.array([0.0, 0.0, 0.0])):
        super(PrimitiveSlater, self).__init__()

        self['n'] = n
        self['l'] = l
        self['m'] = m
        self["exponent"] = exponent
        self["coefficient"] = coefficient
        self['normalization'] = self.normalize()
        self["origin"] = np.array(origin)

    def normalize(self):
        """
        Calculates the normalization constant given by:

        :math:`N_{\lambda p} = [(2n_{\lambda p})!]^{-1/2}(2\zeta_{\lambda p})^{n_{\lambda p}+1/2}`

        """
        output = (
                    (2.0 * self['exponent'])**self['n'] *
                    np.sqrt((2.0 * self['exponent'])/(scipy.misc.factorial(2.0 * self['n'])))
                )

        return output

    def compute(self, coord):
        """
        Calculates the value of PrimitiveSlater at given coordinates (spherical)

        Args:
            coord (array[3]): spherical coordinates (see description below)
            theta (float): [0, 2*pi]; the azimuthal (longitudinal) coordinate.
            phi (float): [0, pi]; the polar (co-latitudinal) coordinate.
            r (float): The r coordinate.

        Returns:
        x (complex float): The value of :math:`\chi_{p \lambda \\alpha}(r, \\theta, \phi)` sampled at ``r``, ``theta`` and ``phi``
        """
        RP = coord - self['origin']

        xy = RP[0]**2 + RP[1]**2
        r = np.sqrt(xy + RP[2]**2)
        theta = np.arctan2(np.sqrt(xy), RP[2])
        phi = np.arctan2(RP[1], RP[0])

        output = (
                self['coefficient'] *
                self['normalization'] *
                self.spherical(theta, phi) *
                self.radial(r)
            )

        return output

    def overlap(self, other):
        """
        Calculates analytically the overlap integral between two SlaterPrimitives.

        Args:
            other (PrimitiveSlater) : function to perform :math:`<\phi_{self} | \phi_{other}>`
        """

        ll = np.abs(self['l'])
        mm = np.abs(self['m'])

        other_ll = np.abs(other['l'])
        other_mm = np.abs(other['m'])

        l = int((float((ll+other_ll+2)-np.abs(ll-other_ll))) / (float((ll+other_ll+2) + np.abs(ll-other_ll))))
        m = int((float((mm+other_mm+2)-np.abs(mm-other_mm))) / (float((mm+other_mm+2) + np.abs(mm-other_mm))))

        n = scipy.misc.factorial(self['n'] + other['n'])

        zeta = self['exponent'] + other['exponent']
        aux = self['n'] + other['n'] + 1

        output = (
                l * m * (n / (zeta**aux)) *
                self['coefficient'] * other['coefficient'] *
                self['normalization'] * other['normalization']
                )

        return output

    def radial(self, r):
        """
        Calculates the radial part of the Slater orbital.

        :math:`R_{\lambda p} = N_{\lambda p} r^{\lambda p -1} e^{-\zeta_{\lambda p} r}`

        where, N is the Normalization given by:

        :math:`N_{\lambda p} = [(2n_{\lambda p})!]^{-1/2}(2\zeta_{\lambda p})^{n_{\lambda p}+1/2}`

        Args:
            r (float): The r coordinate.

        Returns:
            R (float): The radial :math:`R_{\lambda p}` value sampled at ``r``.
        """

        # Radial part eq. 6
        aux = self['n'] - 1
        output = (r**aux) * np.exp(-self['exponent'] * r)

        return output

    def spherical(self, theta, phi):
        """
        Calculates the spherical harmonic part of the Slater orbital. Uses scipy.special.sph_harm

        Args:
            theta (float): [0, 2*pi]; the azimuthal (longitudinal) coordinate.
            phi (float): [0, pi]; the polar (co-latitudinal) coordinate.

        Returns:
            y (float): The harmonic :math:`Y^m_l` sampled at ``theta`` and ``phi``

        """

        # Spherical harmonic eq. 6
        if self['m'] == 0:
            output = sph_harm(self['m'], self['l'], theta, phi)
        elif self['m'] < 0:
            aux_a = sph_harm(self['m'], self['l'], theta, phi)
            aux_b = sph_harm(-self['m'], self['l'], theta, phi)
            output = complex(0.0, 1.0) * 0.70710678118654757 * (aux_b + aux_a)
        else:
            aux_a = sph_harm(self['m'], self['l'], theta, phi)
            aux_b = sph_harm(-self['m'], self['l'], theta, phi)
            output = 0.70710678118654757 * (aux_b - aux_a)

        return output.real

    def show(self):
        """
        Prints information about the object.
        """
        print("    origin: ", self['origin'])
        print("    z: ", self['exponent'])
        print("    coefficient: ", self['coefficient'])
        print("    n: ", self['n'])
        print("    l: ", self['l'])
        print("    m: ", self['m'])
        print("    normalization: ", self['normalization'])
