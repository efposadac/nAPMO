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


def INDEX(i, j):
    if i > j:
        return int((i * (i + 1) / 2) + j)
    else:
        return int((j * (j + 1) / 2) + i)


class PSIN(napmo.PSIB):

    """
    Numerical WaveFunction class

    This class is initiated with some matrices obtained from a minimal basis analytical
    calculation.

    Args:
        psix (PSIX) : Wavefunction object, can be either numerical or analytical
        psi_grid (ndarray) : Wavefunction expanded on the grid
    """

    def __init__(self, psix, grid, ndim=None, debug=False):

        # Initialize base class
        if ndim:
            self._aux_ndim = ndim
        elif psix.species.get('occupation') == 0:
            self._aux_ndim = 1
        else:
            self._aux_ndim = psix.species.get('occupation')

        super(PSIN, self).__init__(psix.species,
                                   ndim=self._aux_ndim)

        self._debug = debug
        self._res = np.ones(self.ndim)
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = psix._pc
        self.species = psix.species
        self._e = psix.O[:self.ndim].copy()
        self._total_mass = psix._total_mass
        self._ecenergy = psix._ecenergy
        self._energy = psix._energy
        self._tf = psix._tf

        if self._tf:
            # TF Correction
            self._mass_inv = ((1.0 / self.species.get('mass')) - (1.0 / self._total_mass))
        else:
            self._mass_inv = 1.0 / self.species.get('mass')

        self._grid = grid
        self._lmax = int(napmo.lebedev_get_order(self._grid.nang) / 2)

        # Calculate the wave-function on the grid for occupied orbitals only
        gbasis = self.species.get('basis').compute(self._grid.points).T.copy()
        self.psi = self._compute_psi_from_cm(psix.C, gbasis)

        # Initialize integrals
        self.Kgrid = np.zeros([self.ndim, self._grid.size])
        self.Jgrid = np.zeros(self._grid.size)
        self.XCgrid = np.zeros(self._grid.size)
        self.Vnuc = napmo.compute_nuclear(self._grid, self._pc)
        self._exchange = False

        # Compute integrals
        self.normalize()
        self.compute_1body()
        self.compute_guess()

        # Initialize optimizer
        self._optimize = napmo.PSIO(self)

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

        # print("\n Overlap Matrix:", self.symbol)
        # print(self.S)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        self._compute_kinetic_operator()
        self.T[:] = self._get_operator_matrix(self.Tgrid)

        # print("\n Kinetic Matrix (N): "+ self.symbol + ": ", self.T.sum())
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """

        self._compute_nuclear_operator()
        self.V[:] = self._get_operator_matrix(self.Vgrid)
        self.V *= self.species.get('charge')

        # print("\n Attraction Matrix", self.symbol)
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix.
        """

        self._compute_density()

        # FELIX: add conditional for HF or hybrid
        if self.species.get('size') > 1:
            self._compute_2body_coulomb()
            if self._exchangefactor != 0.0:
                self._compute_2body_exchange()

            with napmo.runtime.timeblock('Numerical 2 body'):
                napmo.cext.nwavefunction_compute_2body_matrix_mol(
                    byref(self), self._grid._this, self.psi, self.Jgrid, self.Kgrid)

            self.G *= self.species.get('charge')

        # print("\n G Matrix:", self.symbol)
        # print(self.G)

    def compute_coupling(self, other_psi, direct=False):
        """
        Computes the two-body coupling matrix

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
            direct (bool) : Whether to calculate eris on-the-fly or not
        """

        # TODO check Cgrid for more than one species
        aux = np.zeros([self._aux_ndim, self._aux_ndim])
        self.J[:] = 0.0
        for psi in other_psi:
            if self.sid != psi.sid:

                psi.Cgrid = napmo.compute_coulomb(
                    psi._grid, psi.Dgrid.sum(axis=0), psi.lmax)

                psi.Cgrid *= psi.species.get('charge')

                napmo.cext.nwavefunction_compute_coupling(
                    byref(self), self._grid._this, self.psi, psi.Cgrid, aux)

                self.J += aux

        # TODO: Check if J is multiplied many times by the charge when there are three species
        self.J *= self.species.get('charge')

        # print("\n Coupling Matrix: ", self.symbol)
        # print(self.J)

    def compute_exccor(self):
        """
        Computes the exchange correlation matrix

        Args:
        """
        self._ecenergy = 0.0
        self.XC[:] = 0.0
        self.XCgrid[:] = 0.0

        if (self.symbol == "e-"):
            napmo.cext.nwavefunction_compute_exccor_matrix(
            byref(self), self._grid._this, self.psi, self.Dgrid.sum(axis=0), self.XCgrid)

        # print("\n XC Energy:" + self.symbol + ":")
        # print(self._ecenergy)

        # print("\n XC Potential:" + self.symbol + ":")
        # print(self.XCgrid)

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)

    def compute_cor2species(self,other_psi):
        """
        Computes the exchange correlation matrix

        Args:
        """
        for psi in other_psi:
            if self.sid != psi.sid:
                napmo.cext.nwavefunction_compute_cor2species_matrix(
                    byref(self), byref(psi), self._grid._this, self.psi, self.Dgrid.sum(axis=0), psi.Dgrid.sum(axis=0), self.XCgrid)

        # print("\n XC Energy:" + self.symbol + ":")
        # print(self._ecenergy)

        # print("\n XC Matrix:" + self.symbol + ": ")
        # print(self.XC)

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

    def _compute_density(self):
        """
        Compute the density on the grid
        """
        self.Dgrid = self.psi**2 * self._eta

        # Debug information (Suppose to be the # of e-)
        # print('\nDENS on numeric. Number of ' + self.symbol +
        #       ': ', self._grid.integrate(self.Dgrid.sum(axis=0)))

    def _compute_density_from_dm(self, dm, psi):
        """
        Computes the density in the grid for each orbital from density matrix.
        Includes virtual orbitals.
        """
        with napmo.runtime.timeblock('Numerical Density'):
            dens = np.array([phi * dm.dot(phi) for phi in psi.T]).T

        # Debug information (Suppose to be the # of particles)
        # print('\nDENS Init: ' + self.symbol +
        #       ': ', self._grid.integrate(dens.sum(axis=0)))

        return dens

    def _compute_psi_from_cm(self, cm, psi):
        res = np.zeros([self.ndim, self._grid.size])

        with napmo.runtime.timeblock('Numerical PSI'):
            for i in range(self.ndim):
                for j in range(cm.shape[0]):
                    res[i] += cm[j, i] * psi[j]

        # Debug information(Suppose to be the  # of e-)
        # print('\nPSI on numeric. Number of ' + self.symbol + ': ', self._grid.integrate(
        #     res.sum(axis=0) * res.sum(axis=0)) * self._eta)

        return res

    def _compute_kinetic_operator(self):
        """
        Computes the action of the Kinetic operator on a trial function

        :math:`-\dfrac{1}{2} \nabla^{2} \phi_{i}`
        """

        self.Tgrid = np.array([napmo.compute_kinetic(self._grid, phi, self.lmax)
                               for phi in self.psi])

        self.Tgrid *= self._mass_inv

    def _compute_nuclear_operator(self):
        """
        Computes the action of the Nuclear attraction repulsion
        on a trial function

        :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|} \phi_i`

        """

        self.Vgrid = np.array([phi * self.Vnuc for phi in self.psi])

    def _compute_2body_coulomb(self):
        """
        Computes coulomb potential solving Poisson's equation
        following Becke procedure using the current density
        """
        if self.species.get('size') > 0:

            with napmo.runtime.timeblock('Numerical coulomb'):

                self.Jgrid[:] = napmo.compute_coulomb(
                    self._grid, self.Dgrid.sum(axis=0), self.lmax)

                self.Jgrid *= self.species.get('charge')

        # Debug information
        # print("Coulomb energy " + self.symbol + ": ", 0.5 *
        #       self._grid.integrate(self.Jgrid * self.Dgrid.sum(axis=0)))

    def _compute_2body_exchange(self):
        """
        Computes the exchange potential solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2007).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503

        """
        if not self._exchange:
            if self.species.get('size') > 0:

                with napmo.runtime.timeblock('Numerical exchange'):

                    aux = np.array([napmo.compute_coulomb(
                        self._grid, self.psi[s, :] * self.psi[v, :], self.lmax)
                        for s in range(self.ndim)
                        for v in range(self.ndim)
                        if s >= v])

                    self.Kgrid[:] = np.array([np.array([self.psi[j, :] * aux[INDEX(i, j)]
                                                        for j in range(self.ndim)]).sum(axis=0)
                                              for i in range(self.ndim)])

                    self.Kgrid *= self.species.get('charge')

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

    def _compute_potential(self, coupling):
        """
        Computes the total potential
        """
        # print("Vnuc", self.Vnuc.sum() * self.species.get('charge'))
        # print("Jgrid", self.Jgrid.sum() * self.species.get('charge'))
        # print("Kgrid", self.Kgrid.sum() * self.species.get('charge'))
        # print("coupling", coupling.sum() * self.species.get('charge'))
        # print("XCgrid", self.XCgrid.sum())

        return np.array([(self.Vnuc + self.Jgrid + coupling) * self.species.get('charge')+self.XCgrid])

    def _compute_residual(self, coupling):
        """
        Build R (Eq. 10) :math:`R = (T + V - e) \phi`
        """

        V_tot = self._compute_potential(coupling)

        self.Rgrid = np.array([T + ((V - e) * phi) - (K * self.species.get('charge'))
                               for T, V, e, phi, K in zip(self.Tgrid, V_tot, self.O, self.psi, self.Kgrid)])

    def _compute_energy_correction(self):
        """
        Computes energy correction. Eq. 11
        """
        self.delta_e = np.array([self._grid.integrate(r * phi)
                                 for r, phi in zip(self.Rgrid, self.psi)])

        if self._debug:
            print("Delta Energy " + self.symbol + ": ", self.delta_e)

    def _compute_delta_psi(self, coupling):
        """
        Computes \Delta \psi. Eq. 13
        """

        V_tot = self._compute_potential(coupling)

        # print("V_tot", V_tot.sum())

        self.delta_psi = np.array([
            napmo.compute_dpsi(self._grid, self.lmax, phi, doi, oi, ri, V, self._mass_inv, self.species.get('charge'))
            for phi, doi, oi, ri, V in zip(self.psi, self.delta_e, self._e, self.Rgrid, V_tot)])

        # print("psi", self.psi.size)
        # print("delta psi", self.delta_psi.size)

    def _compute_delta_orb(self):
        """
        Computes \Delta_{orb} Eq. 18
        """
        self._res = np.array([self._grid.integrate(psi * psi)**0.5
                              for psi in self.delta_psi])

        if self._debug:
            print("Delta Orbital " + self.symbol + ": ", self._res)

    def optimize_psi(self, scf, other_psi=None):
        """
        Calculates \Delta \phi as Eq. 14 Becke's paper.
        """
        if self._debug:
            print('Optimizing Wavefunction....')

        if np.abs(self._res.sum()) < 1.0e-5:
            print("No wave-function optimization!!!")
        else:
            aux = np.zeros(self.Jgrid.shape)
            if other_psi is not None:
                aux[:] = np.array([psi.Cgrid
                                   for psi in other_psi
                                   if psi.sid != self.sid]).sum(axis=0)

            self._compute_residual(aux)
            self._compute_energy_correction()
            self._compute_delta_psi(aux)
            self._compute_delta_orb()

            # TODO: check optimization
            self.psi, self._e = self._optimize.optimize(self, scf, other_psi)

            self.normalize()
            self.compute_1body()
            self._exchange = False

            if self._debug:
                print('...Done!')

    @property
    def lmax(self):
        return self._lmax
