# file: psi_analytic.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


class PSIA(napmo.PSIB):
    """
    Defines the Fock operator for a Hartree-Fock Calculation with analytic calculation of integrals.
    """

    def __init__(self, species, point_charges):
        super(PSIA, self).__init__(species)

        # Initialize Libint object to calculate integrals
        self._libint = napmo.cext.LibintInterface_new(species.get('id'))
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = point_charges

        for particle in species.get('particles'):
            napmo.cext.LibintInterface_add_basis(
                self._libint, particle.get('basis')._this)

        for point in point_charges:
            napmo.cext.LibintInterface_add_pointcharges(
                self._libint,
                point.get('atomic_number', point.get('charge')),
                point.get('origin'))

        # Initialize all
        self.compute_overlap()
        self.compute_kinetic()
        self.compute_nuclear()
        self.compute_hcore()
        self.compute_guess()

    def compute_overlap(self):
        """
        Computes the overlap matrix
        """
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 1, self.S)
        # print("\n Overlap Matrix:")
        # print(self.S)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 2, self.T)
        self.T /= self.species.get('mass')
        # print("\n Kinetic Matrix:")
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 3, self.V)
        self.V *= -self.species.get('charge')
        # print("\n Attraction Matrix")
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix

        Args:
            direct (bool) : Whether to calculate eris on-the-fly or not
        """
        if self.species.get('size') > 1:

            napmo.cext.LibintInterface_init_2body_ints(self._libint)

            if direct:
                napmo.cext.LibintInterface_compute_2body_direct(
                    self._libint, self.D, self.G)
            else:
                if not self._ints:
                    self._ints = napmo.cext.LibintInterface_compute_2body_ints(
                        self._libint, self.D)
                napmo.cext.wavefunction_compute_2body_matrix(
                    byref(self), self._ints)

            self.G *= self.species.get('charge')**2

            # print("\n G Matrix:")
            # print(self.G)

    def compute_coupling(self, other_psi, direct=False):
        """
        Computes the two-body coupling matrix

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
            direct (bool) : Whether to calculate eris on-the-fly or not
        """
        aux = np.zeros([self.nbasis, self.nbasis])
        self.J[:] = 0.0
        for psi in other_psi:
            if self.sid != psi.sid:
                napmo.cext.LibintInterface_compute_coupling_direct(
                    self._libint, psi._libint, psi.D, aux)
                aux *= self.species.get('charge') * psi.species.get('charge')
                self.J += aux

        self.J += np.tril(self.J, -1).T

        # print("\n Coupling Matrix:")
        # print(self.J)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.H[:] = self.T + self.V
        # print("\n hcore Matrix")
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using the *HCORE* method:
        :math:`H C = e S C`
        """
        napmo.cext.wavefunction_guess_hcore(byref(self))
        # print("\n Coefficients Matrix")
        # print(self.C)

        # print("\n Density Guess")
        # print(self.D)

    def build_fock(self):
        """
        Builds the Fock matrix
        """

        self.F[:] = self.H + self.G + self.J
        # print("\n Fock Matrix:")
        # print(self.F)
