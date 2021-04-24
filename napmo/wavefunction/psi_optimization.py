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

    def __init__(self, psin, psi, ndim=None):
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
        self._psi = psi

    def optimize(self, rho, prev_psi):
        """
        Calculates next psi on the auxiliary basis, it uses the (i-1)th rho (density)
        and orbitals to calculate the new ones.
        """
        self.Dgrid = rho
        self.prev_psi = prev_psi

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
        following the procedure of Shiozaki, T., & Hirata, S. (2007).
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
