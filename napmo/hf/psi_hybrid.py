# file: psi_hybrid.py
# nAPMO package
# Copyright (c) 2016-2017, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo

import sys


class PSIH(napmo.PSIA):
    """
    Defines the Fock operator for a Hartree-Fock Calculation with hybrid calculation of integrals.
    """

    def __init__(self, species, point_charges, total_mass, grid, psia=None):
        super(PSIH, self).__init__(species, point_charges, total_mass)

        self._grid = grid
        self._lmax = int(napmo.lebedev_get_order(self._grid.nang) / 2)
        self.psi = self.species.get('basis').compute(
            self._grid.points).T.copy()

        self._res = []

        if psia:
            self.S[:] = psia.S
            self.T[:] = psia.T
            self.V[:] = psia.V
            self.H[:] = psia.H
            self.C[:] = psia.C
            self.D[:] = psia.D
            self.L[:] = psia.L
            self.G[:] = psia.G
            self.J[:] = psia.J
            self.XC[:] = psia.XC
            self.F[:] = psia.F
            self.O[:] = psia.O

            self._energy = psia._energy
            self._rmsd = psia._rmsd


    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix

        Args:
            direct (bool) : Whether to calculate eris on-the-fly or not
        """

        self._compute_density()

        if self.species.get('size') > 1:
            
            with napmo.runtime.timeblock('Analytical G matrix'):

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

        aux = np.zeros([self.ndim, self.ndim])
        self.J[:] = 0.0

        with napmo.runtime.timeblock('Numerical coupling ints'):
            for psi in other_psi:
                if self.sid != psi.sid:

                    psi.Cgrid = napmo.compute_coulomb(
                        psi._grid, psi.Dgrid.sum(axis=0), psi.lmax)

                    psi.Cgrid *= self.species.get('charge') * \
                        psi.species.get('charge')

                    napmo.cext.nwavefunction_compute_coupling(
                        byref(self), self._grid._this, self.psi, psi.Cgrid, aux)

                    self.J += aux

        # Debug information
        # print("\n Coupling Matrix " + self.symbol + ": ")
        # print(self.J)

    def compute_exccorr(self):
        """
        Computes the exchange correlation matrix

        Args:
        """

        print("\n XC Matrix:" + self.symbol + ": ", self.XC.sum())
        print(self.XC)

        # print("\n D Matrix:" + self.symbol + ": ", self.D.sum())
        # if self.species.get('size') > 0:

        #     with napmo.runtime.timeblock('Numerical coupling ints'):

        #         napmo.cext.LibintInterface_init_2body_ints(self._libint)

        #         if direct:
        #             napmo.cext.LibintInterface_compute_2body_direct(
        #                 self._libint, self.D, self.G)
        #         else:
        #             if not self._ints:
        #                 self._ints = napmo.cext.LibintInterface_compute_2body_ints(
        #                     self._libint, self.D)
        #             napmo.cext.wavefunction_compute_2body_matrix(
        #                 byref(self), self._ints)

        #         self.G *= self.species.get('charge')**2

            # print("\n G Matrix:" + self.symbol + ": ", self.G.sum())
            # print(self.G)


    def _compute_density(self):
        """
        Compute the density on the grid
        """
        with napmo.runtime.timeblock('Numerical Density'):
            self.Dgrid = np.array([phi * self.D.dot(phi)
                                   for phi in self.psi.T]).T

        # Debug information (Suppose to be the # of e-)
        # print('\nDENS on hybrid. Number of ' + self.symbol +
        #       ': ', self._grid.integrate(self.Dgrid.sum(axis=0)))

    def optimize_psi(self, scf, other_psi=None):
        pass


    @property
    def lmax(self):
        return self._lmax
