# file: psi_numeric.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo
import matplotlib.pyplot as plt


class PSIN(napmo.PSIB):

    """
    Numerical WaveFunction class

    This class is initiated with some matrices obtained from a minimal basis analytical
    calculation.

    Args:
        psix (PSIX) : Wavefunction object, can be either numerical or analytical
        psi_grid (ndarray) : Wavefunction expanded on the grid
    """

    def __init__(self, psix, grid):

        # Initialize base class
        super(PSIN, self).__init__(psix.species,
                                   ndim=psix.species.get('occupation'))

        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = psix._pc
        self.species = psix.species
        self._grid = grid

        aux = int(self.ndim * (self.ndim + 1) / 2)
        self.Kgrid = np.zeros([aux, self._grid.size])

        self.Jgrid = np.zeros(self._grid.size)
        self.Vgrid = np.zeros([self.ndim, self._grid.size])

        self._lmax = int(napmo.lebedev_get_order(self._grid.nang) / 2)
        self._e = psix.O[:self.ndim]
        self.Vnuc = napmo.compute_nuclear(self._grid, self._pc)
        self.Vnuc *= self.species.get('charge')

        self._optimize = napmo.PSIO(self)

        # Calculate the wavefunction on the grid for occupied orbitals only
        gbasis = self.species.get('basis').compute(self._grid.points)
        self.psi = self._compute_psi_from_cm(psix.C, gbasis.T)

        # Initialize integrals
        self.normalize()
        self.compute_1body()
        self.compute_guess()
        self._exchange = False

    def normalize(self):
        """
        Normalizes the orbitals on the grid.
        """
        self.compute_overlap()
        norm = 1.0 / np.sqrt(self.S.diagonal())
        self.psi = np.array([psi * n for psi, n in zip(self.psi, norm)])

    def compute_1body(self):
        """
        Computes all 1 body integrals
        """
        with napmo.runtime.timeblock('Numerical 1 body'):
            self.compute_overlap()
            self.compute_kinetic()
            self.compute_nuclear()
            self.compute_hcore()

    def compute_overlap(self):
        """
        Computes the overlap matrix
        """
        self.S[:] = self._get_operator_matrix(self.psi)

        # print("\n Overlap Matrix:")
        # print(self.S)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        self._compute_kinetic_operator()
        self.T[:] = self._get_operator_matrix(self.Tgrid)

        # print("\n Kinetic Matrix:")
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """
        self._compute_nuclear_operator()
        self.V[:] = self._get_operator_matrix(self.Vgrid)

        # print("\n Attraction Matrix")
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix.
        """

        self._compute_density()
        self._compute_2body_coulomb()
        self._compute_2body_exchange()

        with napmo.runtime.timeblock('Numerical 2 body'):
            napmo.cext.nwavefunction_compute_2body_matrix(
                byref(self), self._grid._this, self.psi, self.Jgrid, self.Kgrid)

            self.G *= self.species.get('charge')**2

        # print("\n G Matrix:")
        # print(self.G)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.H[:] = self.T + self.V

        # print("\n hcore Matrix")
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using analytical density.
        """
        napmo.cext.wavefunction_guess_hcore(byref(self))

        # print("\n Guess Matrix")
        # print(self.D)

    def build_fock(self):
        """
        Builds the Fock matrix
        """

        self.F[:] = self.H + self.G + self.J

        # print("\n Fock Matrix:")
        # print(self.F)

    def _compute_density(self):
        """
        Compute the density on the grid
        """

        self.Dgrid = self.psi**2 * self._eta

        # Debug information (Suppose to be the # of e-)
        print('0Number of -e: ', self._grid.integrate(self.Dgrid.sum(axis=0)))

    def _compute_density_from_dm(self, dm, psi):
        """
        Computes the density in the grid for each orbital from density matrix. 
        Includes virtual orbitals.
        """
        with napmo.runtime.timeblock('Numerical Density'):
            dens = np.array([phi * dm.dot(phi)
                             for phi in psi.T]).T

        return dens

    def _compute_psi_from_cm(self, cm, psi):
        res = np.zeros([self.ndim, self._grid.size])

        with napmo.runtime.timeblock('Numerical Density'):
            for i in range(self.ndim):
                for j in range(cm.shape[0]):
                    res[i] += cm[j, i] * psi[j]

        # Debug information (Suppose to be the # of e-)
        print('PSI Number of -e: ', self._grid.integrate(
            res.sum(axis=0) * res.sum(axis=0)) * self._eta)

        return res

    def _compute_kinetic_operator(self):
        """
        Computes the action of the Kinetic operator on a trial function

        :math:`-\dfrac{1}{2} \nabla^{2} \phi_{i}`
        """

        self.Tgrid = np.array([napmo.compute_kinetic(self._grid, phi, self.lmax)
                               for phi in self.psi])

    def _compute_nuclear_operator(self):
        """
        Computes the action of the Nuclear attraction repulsion
        on a trial function

        :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|} \phi_i`

        """

        self.Vgrid[:] = np.array([phi * self.Vnuc for phi in self.psi])

        # if self.Vgrid.sum() != 0.0:
        #     self.Vgrid[:] = (1.0 / 3.0 * Vgrid) + (2.0 / 3.0 * self.Vgrid)
        # else:
        #     self.Vgrid[:] = Vgrid

    def _compute_2body_coulomb(self):
        """
        Computes coulomb potential solving Poisson's equation
        following Becke procedure using the current density
        """
        if self.species.get('size') > 1:

            with napmo.runtime.timeblock('Numerical coulomb'):

                self.Jgrid[:] = napmo.compute_coulomb(
                    self._grid, self.Dgrid.sum(axis=0), self.lmax)

        # Debug information
        # print("Coulomb energy: ", 0.5 *
        #       self._grid.integrate(Jgrid * self.Dgrid.sum(axis=0)))

    def _compute_2body_exchange(self):
        """
        Computes the exchange potential solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2007).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503

        """
        if not self._exchange:
            if self.species.get('size') > 1:

                with napmo.runtime.timeblock('Numerical exchange'):

                    self.Kgrid[:] = np.array([napmo.compute_coulomb(
                        self._grid, self.psi[s, :] * self.psi[v, :], self.lmax)
                        for s in range(self.ndim)
                        for v in range(self.ndim)
                        if s >= v])

                self._exchange = True

    def _get_operator_matrix(self, O):
        """
        Calculates the matrix representation of the operator ``O``
        based on the grid representation of such operator.

        Each element on the matrix corresponds to:

        :math:`\langle \psi_i |O| \psi_j \rangle`

        Args:
            O (ndarray): :math:`O| \psi_j \rangle` calculated on the grid.
        """
        M = np.array([self._grid.integrate(self.psi[i] * O[j])
                      for i in range(self.ndim)
                      for j in range(self.ndim)])

        M = M.reshape([self.ndim, self.ndim])

        return M

    def _compute_residual(self):
        """
        Build R (Eq. 10) :math:`R = (T + V - e) \phi`
        """

        # self.Rgrid = np.array(
        #     [T + (self.Vnuc * psi) + (self.Jgrid * psi) - K - (e * psi)
        # for T, K, e, psi in zip(self.Tgrid, self.Kgrid, self._e, self.psi)])

        # TODO: K missing!
        self.Rgrid = np.array(
            [T + (self.Vnuc * psi) + (0.5 * (self.Jgrid * psi)) - (e * psi)
             for T, e, psi in zip(self.Tgrid, self._e, self.psi)])

    def _compute_energy_correction(self):
        """
        Computes energy correction. Eq. 11
        """

        self.delta_e = np.array([self._grid.integrate(r * phi)
                                 for r, phi in zip(self.Rgrid, self.psi)])

        print("Delta Energy: ", self.delta_e)

    def _compute_delta_psi(self):
        """
        Computes \Delta \psi. Eq. 13
        """

        # TODO: K missing!
        self.delta_psi = np.array([napmo.compute_dpsi(
            self._grid, self.lmax, phi, doi, oi, ri, self.Vnuc, self.Jgrid)
            for phi, doi, oi, ri in zip(self.psi, self.delta_e, self._e,
                                        self.Rgrid)])

    def _compute_delta_orb(self):
        """
        Computes \Delta_{orb} Eq. 18
        """
        self._res = np.array([self._grid.integrate(psi * psi)**0.5
                              for psi in self.delta_psi])

        print("Delta Orbital: ", self._res)

    def optimize_psi(self, scf):
        """
        Calculates \Delta \phi as Eq. 14 Becke's paper.
        """
        self._compute_residual()
        self._compute_energy_correction()
        self._compute_delta_psi()
        self._compute_delta_orb()

        self.psi, self._e = self._optimize.optimize(self, scf)

        self.compute_1body()
        self._exchange = False

    @property
    def lmax(self):
        return self._lmax
