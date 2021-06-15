# file: psi_optimization.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
from numpy import linalg as LA
from scipy.optimize import minimize

from ctypes import *

import napmo


class PSIO(napmo.PSIN):
    """
    Class for wavefunction; Calculates orbitals and density on the auxiliary basis
    """

    def __init__(self, psin, psi=None, ndim=None):
        if ndim is not None:
            aux = ndim
        else:
            aux = psin.ndim

        super(napmo.PSIN, self).__init__(psin.species, ndim=aux)

        self._debug = psin._debug
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = psin._pc
        self._grid = psin._grid
        self._lmax = psin._lmax
        self._mass_inv = psin._mass_inv
        self._exchange = False

        self._Vnuc = psin.Vnuc
        self._iterations = 1

        # psi here is the auxiliary basis.
        if psi is not None:
            self._psi = psi

        self.Cgrid = 0.0

    def compute_coupling_operator(self, other_psi):
        """
        Computes the two-body coupling matrix

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
        """
        aux = np.array([
            psi.Jpot
            for psi in other_psi
            if psi.sid != self.sid
        ]).sum(axis=0)

        self.Cgrid = np.array([
            aux * self.psi[k]
            for k in range(self.ndim)
        ])

    def build_fock(self):
        """
        Builds the Fock matrix
        """
        self.F[:] = self.get_operator_matrix(self.Fgrid)

        # print("\n Fock Matrix:", self.symbol)
        # print(self.F)

    def compute_xc_grid(self):
        pass
        """
        Computes the exchange correlation matrix

        Args:
        """
        # self._ecenergy = 0.0
        # self.XC[:] = 0.0
        # self.XCgrid[:] = 0.0

        # if (self.symbol == "e-"):
        #     napmo.cext.nwavefunction_compute_exccor_matrix(
        #         byref(self), self._grid._this, self.psi, self.Dgrid.sum(axis=0), self.XCgrid)

        # # print("\n XC Energy:" + self.symbol + ":")
        # # print(self._ecenergy)

        # # print("\n XC Potential:" + self.symbol + ":")
        # # print(self.XCgrid)

        # # print("\n XC Matrix:" + self.symbol + ": ")
        # # print(self.XC)

    def compute_c_2species_grid(self, other_psi):
        pass
        """
        Computes the exchange correlation matrix

        Args:
        """
        # for psi in other_psi:
        #     if self.sid != psi.sid:
        #         napmo.cext.nwavefunction_compute_cor2species_matrix(
        #             byref(self), byref(psi), self._grid._this, self.psi, self.Dgrid.sum(axis=0), psi.Dgrid.sum(axis=0), self.XCgrid)

        # print("\n XC Energy:" + self.symbol + ":")
        # print(self._ecenergy)

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)

    @property
    def iterations(self):
        return self._iterations


class PSIO_BECKE(PSIO):
    """
    Wavefunction optimization using a Becke-based scheme.
    See: Becke, A. D., & Dickson, R. M. (1990).
    Numerical solution of Schrödinger’s equation in polyatomic molecules.
    The Journal of Chemical Physics, 92(6), 3610–3612. https://doi.org/10.1063/1.457869
    """

    def __init__(self, psin):
        super(PSIO_BECKE, self).__init__(psin, ndim=psin.ndim * 2)
        self.prev_O = psin._prev_O
        self.delta_orb = np.zeros(self.ndim)
        self.scf = napmo.SCF({"print": True})

    def optimize(self, psi, other_psi=None):
        self.compute_potential(psi)
        self.compute_residual(psi)
        self.compute_energy_correction(psi)
        self.compute_delta_psi(psi)

        # Variational opt
        self._psi = np.vstack([psi.psi, self.delta_psi])
        self.Dgrid = psi.Dgrid

        if other_psi is not None:
            self.compute_coupling_operator(other_psi)

        # Compute Fock
        self.compute_overlap()
        self.compute_1body_operator()
        self.compute_coulomb_operator(update_rho=False)
        self.compute_exchange_operator()
        self.build_fock_operator()
        self.build_fock()

        # Symmetrizing Fock
        self.asymmetry = (np.abs(self.F - self.F.T)).max()

        if self.iterations == 1:
            print("Asymmetry in Fock Matrix:", self.asymmetry)
            print("Warning: Fock Matrix symmetrized!")

        self.F[:] = (self.F + self.F.T) / 2.0

        napmo.cext.wavefunction_iterate(byref(self))
        # self.scf.single(self, pprint=True)
        self.compute_delta_orb(psi)

        # New Orbitals
        self.prev_O = psi.O
        self._iterations += 1
        self.next_psi = self.compute_psi_from_cm(self.C, self.psi)

        return self.next_psi

    def compute_potential(self, psi):
        self.V_tot = np.array([
            v + (p * psi.Jpot) + (k * psi._x_factor)
            for t, v, p, k in zip(psi.Tgrid, psi.Vgrid, psi.psi, psi.Kpot)
        ])

    def compute_residual(self, psi):

        self.Rgrid = np.array([
            T + V - (e * phi)
            for T, V, e, phi in zip(psi.Tgrid, self.V_tot, self.prev_O, psi.psi)
        ])

    def compute_energy_correction(self, psi):
        self.delta_e = np.array([
            psi.grid.integrate(r * phi)
            for r, phi in zip(self.Rgrid, psi.psi)
        ])

    def compute_delta_psi(self, psi):
        self.delta_psi = np.array([
            napmo.compute_dpsi(psi.grid, psi.lmax, phi, doi, oi, ri, V, psi._mass_inv)
            for phi, doi, oi, ri, V in zip(psi.psi, self.delta_e, self.prev_O, self.Rgrid, self.V_tot)
        ])

    def compute_delta_orb(self, psi):
        self.delta_orb = self.O[:psi.ndim] - psi.O
        # print(np.array([
        #     np.sqrt(psi.grid.integrate(p * p))
        #     for p in self.delta_psi
        # ]))


class PSIO_SO(PSIO):
    """
    Wavefunction optimization using an auxiliary basis proposed by So Hirata.
    """

    def __init__(self, psin):
        self._aobasis = napmo.AuxiliaryBasis(psin.grid)
        self.delta_orb = np.zeros(self.ndim)
        super(PSIO_SO, self).__init__(psin, psi=self.aobasis.basis, ndim=self.aobasis.nao)

    def optimize(self, psi, other_psi=None):
        """
        Calculates next psi on the auxiliary basis, it uses the (i-1)th rho (density)
        and orbitals to calculate the new ones.
        """
        self.Dgrid = psi.Dgrid
        self.prev_psi = psi.psi

        if other_psi is not None:
            self.compute_coupling_operator(other_psi)

        # Compute Fock
        if self.iterations == 1:
            self.compute_overlap()
            self.compute_1body_operator()

        self.compute_coulomb_operator(update_rho=False)
        self.compute_exchange_operator()
        self.build_fock_operator()
        self.build_fock()

        # Symmetrizing Fock
        self.asymmetry = (np.abs(self.F - self.F.T)).max()

        if self.iterations == 1:
            print("Asymmetry in Fock Matrix:", self.asymmetry)
            print("Warning: Fock Matrix symmetrized!")

        self.F[:] = (self.F + self.F.T) / 2.0

        napmo.cext.wavefunction_compute_coefficients(byref(self))

        # New Orbitals
        self.next_psi = self.compute_psi_from_cm(self.C, self.psi)

        self._iterations += 1

        return self.next_psi

    def compute_exchange_operator(self):
        """
        Computes the exchange potential solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2017).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503
        """
        with napmo.runtime.timeblock('Numerical exchange opt'):
            self.Kgrid = np.array([
                np.array([
                    self.prev_psi[s] * napmo.compute_coulomb(
                        self.grid, self.prev_psi[s] * self.psi[v], self.lmax
                    )
                    for s in range(self.occupation)
                ]).sum(axis=0)
                for v in range(self.ndim)
            ]) * self.species.get('charge')

    @property
    def aobasis(self):
        return self._aobasis
