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
from scipy import linalg as SLA


class PSIA(napmo.PSIB):
    """
    Defines the Fock operator for a Hartree-Fock Calculation with analytic calculation of integrals.
    """

    def __init__(self, species, point_charges, total_mass, grid=None, options=None):
        super(PSIA, self).__init__(species, options=options)

        # Initialize Libint object to calculate integrals
        self._libint = napmo.cext.LibintInterface_new(species.get('id'))
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = point_charges
        self._total_mass = total_mass
        self._grid = grid

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
        self.compute_transformation()
        self.compute_kinetic()
        self.compute_nuclear()
        self.compute_hcore()
        self.compute_guess()

        if grid is not None:
            self._gbasis = self.species.get('basis').compute(self._grid.points).T.copy()

    def compute_overlap(self):
        """
        Computes the overlap matrix
        """
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 1, self.S)
        # print("\n Overlap Matrix:" + self.symbol + ": ", self.S.sum())
        # print(self.S)

    def compute_transformation(self):
        """
        Computes the overlap matrix
        """
        # S = np.fromfile("Overlap.txt", sep=" ", dtype=np.float64)
        # S = S.reshape(self.ndim, self.ndim)
        # self.S[:] = S

        tmp = SLA.lapack.dsyev(self.S)
        self.O[:] = tmp[0]
        self.X[:] = tmp[1]
        napmo.cext.wavefunction_transformation_matrix(byref(self))
        # print("\n Transformation Matrix:" + self.symbol + ": ", self.X.sum())
        # print(self.X)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 2, self.T)

        if self._tf:
            # Translational correction
            self.T *= ((1.0 / self.species.get('mass')) - (1.0 / self._total_mass))
        else:
            self.T /= self.species.get('mass')

        # print("\n Kinetic Matrix (A): "+ self.symbol + ": ", self.T.sum())
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """
        if len(self._pc) > 0:
            napmo.cext.LibintInterface_compute_1body_ints(
                self._libint, 3, self.V)

        self.V *= -self.species.get('charge')
        # print("\n Attraction Matrix" + self.symbol + ": ", self.V.sum())
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix

        Args:
            direct (bool) : Whether to calculate eris on-the-fly or not
        """

        # print("\n D Matrix:" + self.symbol + ": ", self.D.sum())

        if self.species.get('size') > 0:
            with napmo.runtime.timeblock('Numerical coupling ints'):

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

            # print("\n G Matrix:" + self.symbol + ": ", self.G.sum())
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

        # print("\n Coupling Matrix " + self.symbol + ": ", self.J.sum())
        # print(self.J)

    def compute_xc(self):
        """
        Computes the exchange correlation matrix - numerically
        """
        if self._functional is not None:
            # Density in the grid
            self.Dgrid = self._compute_density_from_dm(self.D, self.gbasis)
            cene, cpot = func.compute_correlation(self.Dgrid)
            xene, xpot = func.compute_exchange(self.Dgrid)

            print(cene, cpot)

            # self._xc_energy = 0.0
            # self.XC[:] = 0.0


            # # Exchange correlation potential in the grid - Probably it's not required in the analytical SCF
            # XCgrid = np.zeros(grid.size)

            # if (self.symbol == "e-"):
            #     napmo.cext.nwavefunction_compute_exccor_matrix(
            #         byref(self), grid._this, gbasis, self.Dgrid.sum(axis=0), XCgrid)

        # print("\n XC Energy:" + self.symbol + ":")
        # print(self._ecenergy)

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)

    # def compute_cor2species(self, other_psi):
    #     """
    #     Computes the exchange correlation matrix - numerically

    #     Args:
    #         other_psi (WaveFunction) : WaveFunction object for the other species.
    #     """

    #     # numerical wavefunction
    #     # FELIX: TODO, use input grid
    #     grid = napmo.BeckeGrid(self.species, 100, 110)

    #     # Atomic orbitals represented in the grid
    #     gbasis = self.species.get('basis').compute(grid.points).T.copy()

    #     # Density in the grid
    #     self.Dgrid = self._compute_density_from_dm(self.D, gbasis)

    #     # Exchange correlation potential in the grid - Probably it's not required in the analytical SCF
    #     XCgrid = np.zeros(grid.size)

    #     for psi in other_psi:
    #         if self.sid != psi.sid:
    #             othergbasis = psi.species.get('basis').compute(grid.points).T.copy()
    #             psi.Dgrid = psi._compute_density_from_dm(psi.D, othergbasis)
    #             napmo.cext.nwavefunction_compute_cor2species_matrix(
    #                 byref(self), byref(psi), grid._this, gbasis, self.Dgrid.sum(axis=0), psi.Dgrid.sum(axis=0), XCgrid)

        # print("\n XC Energy:" + self.symbol + ":")
        # print(self._ecenergy)

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.H[:] = self.T + self.V
        # print("\n hcore Matrix" + self.symbol + ": " , self.H.sum())
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using the *HCORE* method:
        :math:`H C = e S C`
        """
        napmo.cext.wavefunction_guess_hcore(byref(self))
        # print("\n Coefficients Matrix" + self.symbol + ": ", self.C.sum())
        # print(self.C)

        # print("\n Density Guess" + self.symbol + ": ", self.D.sum())
        # print(self.D)

    def build_fock(self):
        """
        Builds the Fock matrix
        """

        self.F[:] = self.H + self.G + self.J + self.XC
        # print("\n Fock Matrix:" + self.symbol + ": ", self.F.sum())
        # print(self.F)
