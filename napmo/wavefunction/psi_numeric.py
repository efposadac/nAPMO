# file: psi_numeric.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


def INDEX(i, j):
    if i > j:
        return int((i * (i + 1) / 2) + j)
    else:
        return int((j * (j + 1) / 2) + i)


class PSIN(napmo.PSIB):

    """
    Numerical WaveFunction class

    This class is initiated with some matrices obtained from a minimal basis
    analytical calculation.

    Args:
        psix (PSIX) : Wavefunction object, can be either numerical or analytical
        psi_grid (ndarray) : Wavefunction expanded on the grid
    """

    def __init__(self, psix, grid, ndim=None, debug=False, aux_basis=False):

        # Initialize base class
        if ndim:
            self._ndim = ndim
        elif psix.species.get('occupation') == 0:
            self._ndim = 1
        else:
            self._ndim = psix.species.get('occupation')

        super(PSIN, self).__init__(psix.species, ndim=self.ndim)

        # Initialization
        self._debug = debug
        self._pc = psix._pc
        self._total_mass = psix._total_mass
        self._diis = napmo.cext.LibintInterface_diis_new(2)

        # TF Correction
        if self._tf:
            self._mass_inv = ((1.0 / self.species.get('mass')
                               ) - (1.0 / self._total_mass))
        else:
            self._mass_inv = 1.0 / self.species.get('mass')

        self._xc_energy = psix._xc_energy
        self._energy = psix._energy
        self._prev_O = psix.O

        self._grid = grid
        self._lmax = int(napmo.lebedev_get_order(self.grid.nang) / 2)

        # Calculate the point charges potential
        self._Vnuc = napmo.compute_nuclear(self.grid, self._pc)

        # Calculate the wave-function on the grid for occupied orbitals only
        self._gbasis = self.species.get(
            'basis').compute(self.grid.points).T.copy()

        self._psi = self.compute_psi_from_cm(psix.C, self._gbasis)

        # Start WF calculation
        self.initialize()
        self.Cgrid = np.zeros(self.grid.size)
        self.Jpot = np.zeros(self.grid.size)
        self.Kpot = np.zeros(self.grid.size)

        # Calculate two-body
        self.calculate_2body = self.species.get('size') > 0 or (
            self.species.get('is_electron') and self.species.get('size') > 0)

        # Calculate Aux Basis
        if aux_basis:
            self._optimize = napmo.PSIO_SO(self)
        else:
            self._optimize = napmo.PSIO_BECKE(self)


    def initialize(self):
        """
        Initialize integrals
        """
        self.normalize()
        self.compute_1body()
        self.compute_guess()

    def normalize(self):
        """
        Normalizes the orbitals on the grid.
        """
        self.compute_overlap()
        norm = 1.0 / np.sqrt(self.S.diagonal())
        self._psi = np.array([psi * n for psi, n in zip(self.psi, norm)])

        # Overlap normalized
        self.compute_overlap()

    def compute_1body(self):
        """
        Computes all 1 body integrals
        """
        with napmo.runtime.timeblock('Numerical 1 body'):
            self.compute_kinetic()
            self.compute_nuclear()
            self.compute_hcore()

    def compute_overlap(self):
        """
        Computes the overlap matrix
        """
        self.S[:] = self.get_operator_matrix(self.psi)

        # print("\n Overlap Matrix:", self.symbol)
        # print(self.S)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        self.compute_kinetic_operator()
        self.T[:] = self.get_operator_matrix(self.Tgrid)

        # print("\n Kinetic Matrix (N): "+ self.symbol + ": ", self.T.sum())
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """

        self.compute_nuclear_operator()
        self.V[:] = self.get_operator_matrix(self.Vgrid)

        # print("\n Attraction Matrix", self.symbol)
        # print(self.V)

    def compute_2body(self, direct=False, update_rho=True):
        """
        Computes the two-body matrix.
        """

        # FELIX: add conditional for HF or hybrid
        if self.calculate_2body:

            self.compute_coulomb_potential(update_rho=update_rho)

            if self._x_factor != 0.0:
                self.compute_exchange_potential()

            with napmo.runtime.timeblock('Numerical G matrix'):
                napmo.cext.nwavefunction_compute_2body_matrix(
                    byref(self), self.grid._this, self.psi, self.Jpot, self.Kpot)

            self.G *= self.species.get('charge')**2

        # print("\n G Matrix:", self.symbol)
        # print(self.G)

    def compute_coulomb_potential(self, update_rho=True):
        """
        Computes the two-body coulomb potential in the grid.
        """
        if update_rho:
            self.compute_density()

        with napmo.runtime.timeblock('Numerical coulomb'):
            self.Jpot = napmo.compute_coulomb(
                self.grid, self.Dgrid.sum(axis=0), self.lmax
            )

    def compute_exchange_potential(self):
        """
        Computes the exchange potential solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2017).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503

        Note:
            This is the contribution to build the G matrix, for the grid operator
            use compute_exchange_operator instead

        """
        with napmo.runtime.timeblock('Numerical exchange'):
            aux = np.array([
                napmo.compute_coulomb(
                    self.grid, self.psi[s] * self.psi[v], self.lmax
                )
                for s in range(self.ndim)
                for v in range(self.ndim)
                if s >= v])

            self.Kpot = np.array([
                np.array([
                    self.psi[k] * aux[int(INDEX(l, j))] * self.D[k, l]
                    for k in range(self.ndim)
                    for l in range(self.ndim)
                ]).sum(axis=0)
                for j in range(self.ndim)
            ])

    def compute_coupling(self, other_psi, direct=False):
        """
        Computes the two-body coupling matrix

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
            direct (bool) : Whether to calculate eris on-the-fly or not
        """

        # TODO: check Cgrid for more than one species
        aux = np.zeros([self.ndim, self.ndim])
        self.J[:] = 0.0
        for psi in other_psi:
            if self.sid != psi.sid:

                napmo.cext.nwavefunction_compute_coupling(
                    byref(self), self.grid._this, self.psi, psi.Jpot, aux)

                self.J += (aux * psi.species.get('charge'))

        # TODO: Check if J is multiplied many times by the charge when there are three species
        self.J *= self.species.get('charge')

        # print("\n Coupling Matrix: ", self.symbol)
        # print(self.J)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.H[:] = self.T + self.V

        # print("\n hcore Matrix", self.symbol)
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using analytical density.
        """
        napmo.cext.wavefunction_guess_hcore(byref(self))
        self.Dgrid = self._compute_density_from_dm(self.D, self.psi)

        # print("\n Guess Matrix", self.symbol)
        # print(self.D)

    def build_fock(self):
        """
        Builds the Fock matrix
        """
        self.F[:] = self.H + self.G + self.J + self.XC

        # print("\n Fock Matrix:", self.symbol)
        # print(self.F)

    def compute_1body_operator(self):
        """
        Computes all 1 body integrals
        """
        with napmo.runtime.timeblock('Numerical 1 body'):
            self.compute_kinetic_operator()
            self.compute_nuclear_operator()
            self.compute_hcore_operator()

    def compute_kinetic_operator(self):
        """
        Computes the action of the Kinetic operator on a trial function

        :math:`-\\dfrac{1}{2} \nabla^{2} \\phi_{i}`
        """
        with napmo.runtime.timeblock('Numerical Kinetic'):
            self.Tgrid = np.array([napmo.compute_kinetic(self.grid, phi, self.lmax)
                                   for phi in self.psi[:self.ndim]])

            self.Tgrid *= self._mass_inv

    def compute_nuclear_operator(self):
        """
        Computes the action of the Nuclear attraction repulsion
        on a trial function

        :math:`\\sum_{A} \\dfrac{Z_{A}}{|r_i -R_A|} \\phi_i`

        """
        with napmo.runtime.timeblock('Numerical Nuclear'):
            self.Vgrid = np.array(
                [phi * self.Vnuc for phi in self.psi[:self.ndim]]
            ) * self.species.get('charge')

    def compute_coulomb_operator(self, update_rho=True):
        """
        Computes coulomb potential solving Poisson's equation
        following Becke procedure using the current density
        """

        if self.calculate_2body:
            self.compute_coulomb_potential(update_rho=update_rho)

            with napmo.runtime.timeblock('Numerical coulomb'):

                self.Jgrid = np.array([
                    psi * self.Jpot for psi in self.psi[:self.ndim]]
                )

        # Debug information
        # print("Coulomb energy " + self.symbol + ": ", 0.5 *
        #       self.grid.integrate(self.Jgrid * self.Dgrid.sum(axis=0)))

    def compute_exchange_operator(self):
        """
        Computes the exchange potential on the grid solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2017).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503

        """
        if self.calculate_2body:
            with napmo.runtime.timeblock('Numerical exchange'):
                self.Kgrid = np.array([
                    np.array([
                        self.psi[s] * napmo.compute_coulomb(
                            self.grid, self.psi[s] * self.psi[v], self.lmax
                        )
                        for s in range(self.ndim)
                    ]).sum(axis=0)
                    for v in range(self.ndim)
                ])

    def compute_hcore_operator(self):
        """
        Builds the Hcore matrix
        """
        self.Hgrid = self.Tgrid + self.Vgrid

        # print("\n hcore Matrix", self.symbol)
        # print(self.H)

    def build_fock_operator(self):
        """
        Builds the Fock matrix expanded on the grid.
        """

        self.Fgrid = self.Hgrid + self.Jgrid - self.Kgrid + self.Cgrid

    def compute_density(self):
        """
        Compute the density on the grid
        """
        self.Dgrid = self.psi * self.psi * self._eta

        # Debug information (Suppose to be the # of e-)
        # print('\nDENS on numeric. Number of ' + self.symbol
        #     + ': ', self.grid.integrate(self.Dgrid.sum(axis=0)))

    def compute_density_from_dm(self, dm, psi):
        """
        Computes the density in the grid for each orbital from density matrix.
        Includes virtual orbitals.
        """
        with napmo.runtime.timeblock('Numerical Density'):
            dens = np.array([phi * dm.dot(phi) for phi in psi.T]).T

        # Debug information (Suppose to be the # of particles)
        # print('\nDENS Init: ' + self.symbol +
        #       ': ', self.grid.integrate(dens.sum(axis=0)))

        return dens

    def compute_psi_from_cm(self, cm, psi):

        with napmo.runtime.timeblock('Numerical PSI'):
            res = np.array([
                (psi.T * cm[:, j]).T.sum(axis=0) for j in range(cm.shape[0])
            ])

        # Debug information(Suppose to be the  # of e-)
        # print('\nPSI on numeric. Number of ' + self.symbol + ': ', self.grid.integrate(
        #     res[:self.occupation].sum(axis=0)**2) * self._eta)

        return res

    def optimize_psi(self, other_psi=None):
        """
        Uses the auxiliary basis to calculate an optimized version of the orbitals.
        """
        self._psi[:] = self.optimize.optimize(self, other_psi=other_psi)[:self.ndim]
        self.normalize()
        self.compute_1body()

    def get_operator_matrix(self, Operator):
        """
        Calculates the matrix representation of the operator ``O``
        based on the grid representation of such operator.

        Each element on the matrix corresponds to:

        :math:`\\langle \\psi_i |O| \\psi_j \rangle`

        Args:
            Operator (ndarray): :math:`O| \\psi_j \rangle` calculated on the grid.
        """
        M = np.array([self.grid.integrate(self.psi[i] * Operator[j])
                      for i in range(self.ndim)
                      for j in range(self.ndim)])

        M = M.reshape([self.ndim, self.ndim])

        return M

    @property
    def grid(self):
        return self._grid

    @property
    def psi(self):
        return self._psi

    @property
    def Vnuc(self):
        return self._Vnuc

    @property
    def lmax(self):
        return self._lmax

    @property
    def optimize(self):
        return self._optimize

    def compute_xc_grid(self):
        #     """
        #     Computes the exchange correlation matrix

        #     Args:
        #     """
        #     self._xc_energy = 0.0
        #     self.XC[:] = 0.0
        #     self.XCgrid[:] = 0.0

        #     if (self.symbol == "e-"):
        #         napmo.cext.nwavefunction_compute_exccor_matrix(
        #         byref(self), self.grid._this, self.psi, self.Dgrid.sum(axis=0), self.XCgrid)

        #     # print("\n XC Energy:" + self.symbol + ":")
        #     # print(self._xc_energy)

        #     # print("\n XC Potential:" + self.symbol + ":")
        #     # print(self.XCgrid)

        #     # print("\n XC Matrix:" + self.symbol + ": ")
        #     # print(self.XC)
        pass

    def compute_c_2species_grid(self, other_psi):
        #     """
        #     Computes the exchange correlation matrix

        #     Args:
        #     """
        #     for psi in other_psi:
        #         if self.sid != psi.sid:
        #             napmo.cext.nwavefunction_compute_cor2species_matrix(
        #                 byref(self), byref(psi), self.grid._this, self.psi, self.Dgrid.sum(axis=0),
        #                 psi.Dgrid.sum(axis=0), self.XCgrid
        #             )

        #     # print("\n XC Energy:" + self.symbol + ":")
        #     # print(self._xc_energy)

        #     # print("\n XC Matrix:" + self.symbol + ": ")
        #     # print(self.XC)
        pass

    def compute_xc_matrix(self):
        """
        Builds the exchange-correlation matrix
        """

        if self._xc_vrho is not None:
            napmo.cext.nwavefunction_compute_xc_matrix(
                byref(self),
                self.grid._this,
                self._gbasis,
                self._xc_vrho
            )

            # reset the grid
            self._xc_vrho[:] = 0.0

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)
