# file: psi_analytic.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo
from scipy import linalg as SLA


def print_matrix(A):
    n = A.shape[0]
    for i in range(n):
        for j in range(n):
            print("%12.6f" % (A[i, j]), end="")
        print("")
    print("")


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
            self._xc_vrho = np.zeros(self._grid.size)

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

        # print("\n Kinetic Matrix (A): "+ self.symbol + ": ", self.T.sum(), self.species.get('mass'))
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """
        if len(self._pc) > 0:
            napmo.cext.LibintInterface_compute_1body_ints(
                self._libint, 3, self.V)

        self.V *= -self.species.get('charge')

        # print("\n Attraction Matrix " + self.symbol + ": ", self.V.sum(), self.species.get('charge'))
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix

        Args:
            direct (bool) : Whether to calculate eris on-the-fly or not
        """

        if self.species.get('size') > 1 or self.species.get('is_electron'):

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

    def compute_xc_grid(self, beta_psi=None):
        """
        Computes the exchange correlation - numerically only electrons
        """

        if self._symbol != "e-beta":
            self._xc_energy = 0.0

        if self._functional is not None:

            # Density in the grid
            ndim = 2 if beta_psi else 1
            rho = np.zeros((ndim * self._grid.size))

            # Alpha set
            self.Dgrid = np.array([phi * self.D.dot(phi) for phi in self._gbasis.T]).T
            rho[::ndim] = self.Dgrid.sum(axis=0)

            # Beta set
            if beta_psi is not None:
                beta_psi.Dgrid = np.array([phi * beta_psi.D.dot(phi) for phi in beta_psi._gbasis.T]).T
                rho[1::ndim] = beta_psi.Dgrid.sum(axis=0)

            # Compute functionals
            c_zk, c_vrho = self._functional.compute_correlation(rho)
            x_zk, x_vrho = self._functional.compute_exchange(rho)

            # Calculate XC Energy
            self._xc_energy = self._grid.integrate(rho[::ndim] * (c_zk[0, :] + x_zk[0, :]))

            if beta_psi is not None:
                beta_psi._xc_energy = beta_psi._grid.integrate(rho[1::ndim] * (c_zk[0, :] + x_zk[0, :]))

            # Save potentials

            # Alpha set
            self._xc_vrho[:] = c_vrho[0, :] + x_vrho[0, :]

            # Beta set
            if beta_psi is not None:
                beta_psi._xc_vrho[:] = c_vrho[1, :] + x_vrho[1, :]

            # print("\n XC Energy: " + self.symbol + ":", self._xc_energy)
            # if beta_psi is not None:
            #     print("\n XC Energy: " + beta_psi.symbol + ":", beta_psi._xc_energy)

    def compute_c_2species_grid(self, other_psi):
        """
        Computes the exchange correlation matrix - numerically

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
        """

        for psi in other_psi:
            if self.sid < psi.sid:
                functional = self.options.get('functional', None)

                if functional is not None:

                    # Get the functional name
                    key = ":".join([self._symbol, psi._symbol])
                    functional = functional.get(key,
                                                self.options.get('functional').get(":".join(key.split(":")[::-1])))

                    # Get the functional function
                    functional = napmo.isc_functional_selector(functional)

                    if functional is not None:

                        # Get common points index
                        index = self.options.get('grids').get_common_points(self.symbol, psi.symbol)

                        # Calculate densities
                        self.Dgrid = np.array([phi * self.D.dot(phi) for phi in self._gbasis.T]).T
                        psi.Dgrid = np.array([phi * psi.D.dot(phi) for phi in psi._gbasis.T]).T

                        # Build rho
                        rho = np.take(self.Dgrid.sum(axis=0), index[:, 0], axis=0)
                        other_rho = np.take(psi.Dgrid.sum(axis=0), index[:, 1], axis=0)

                        # Compute functional
                        c_zk, c_vrho, c_other_vrho = functional(rho, other_rho)

                        # Calculate XC Energy (Can't use grid.integrate)
                        aux = (psi._grid.becke_weights[index[:, 1]] * psi._grid.weights[index[:, 1]] * c_zk * rho).sum()

                        self._xc_energy += aux
                        psi._xc_energy += aux

                        # Update the potentials on grids
                        self._xc_vrho[index[:, 0]] += c_vrho
                        psi._xc_vrho[index[:, 1]] += c_other_vrho

                        # print("\n XC Energy: " + self.symbol + "/" + psi.symbol + ": ", aux)

                    functional = None

    def compute_xc_matrix(self):
        """
        Builds the exchange-correlation matrix
        """

        if self._xc_vrho is not None:
            napmo.cext.nwavefunction_compute_xc_matrix(
                byref(self),
                self._grid._this,
                self._gbasis,
                self._xc_vrho
            )

            # reset the grid
            self._xc_vrho[:] = 0.0

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print_matrix(self.XC)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.H[:] = self.T + self.V

        # print("\n hcore Matrix" + self.symbol + ": ", self.H.sum())
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using the *HCORE* method:
        :math:`H C = e S C`
        """
        napmo.cext.wavefunction_guess_hcore(byref(self))

        # print("\n Coefficients Matrix" + self.symbol + ": ", self.C.sum())
        # print(self.C)

        # print("\n Density Guess " + self.symbol + ": ", self.D.sum())
        # print(self.D)

    def build_fock(self):
        """
        Builds the Fock matrix
        """

        self.F[:] = self.H + self.G + self.J + self.XC

        # print("\n Fock Matrix:" + self.symbol + ": ", self.F.sum())
        # print(self.F)
