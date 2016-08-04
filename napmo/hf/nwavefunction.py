# file: nwavefunction.py
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


def print_matrix(m):
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            print('{0:>7.4f} '.format(m[i, j]), end=' ')
        print('')
    print('')


def test_integration(grid, nbasis, psi, op, an):
    test = np.zeros([nbasis, nbasis])
    for i in range(nbasis):
        for j in range(nbasis):
            test[i, j] = grid.integrate(psi[i, :] * op[j])

    print('\nNumerical:')
    print_matrix(test)

    print('\nAnalytical:')
    print_matrix(an)

    print(test.sum(axis=1))


def index2(i, j):
    if i > j:
        ij = i * (i + 1) / 2 + j
    else:
        ij = j * (j + 1) / 2 + i

    return int(ij)


class NWaveFunction(Structure):

    """
    Numerical WaveFunction class

    This class is initiated with some matrices obtained from a minimal basis analytical
    calculation.

    Args:
        psi (WaveFunction) : Optimized wavefunction object
    """

    _fields_ = [
        ("_S", POINTER(c_double)),  # Overlap
        ("_T", POINTER(c_double)),  # Kinetic
        ("_V", POINTER(c_double)),  # Nuclear
        ("_H", POINTER(c_double)),  # Hcore
        ("_C", POINTER(c_double)),  # Coefficients
        ("_D", POINTER(c_double)),  # Density
        ("_L", POINTER(c_double)),  # Last Density
        ("_G", POINTER(c_double)),  # 2 Body
        ("_J", POINTER(c_double)),  # Coupling
        ("_F", POINTER(c_double)),  # Fock
        ("_O", POINTER(c_double)),  # Orbitals
        ("_nbasis", c_int),
        ("_occupation", c_int),
        ("_eta", c_double),
        ("_kappa", c_double),
        ("_energy", c_double),
        ("_rmsd", c_double)  # Root medium square deviation, for D matrix
    ]

    def __init__(self, psi, grid):
        super(NWaveFunction, self).__init__()

        # Initialization
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = psi._pc
        self._occupation = psi._occupation
        self._eta = psi._eta
        self._kappa = psi._kappa
        self._symbol = psi.symbol
        self._sid = psi.sid
        self._energy = 0.0
        self._rmsd = 1.0
        self._pce = 0.0
        self.species = psi.species
        self._grid = grid
        self._lmax = int(napmo.lebedev_get_order(self._grid.nang) / 2)

        # Calculating initial numerical wavefunction (From basis)
        aux = self.species.get('basis').compute(self._grid.points).T
        self.psi = aux.copy()
        self._nbasis = self.psi.shape[0]

        self.S = np.zeros([self.nbasis, self.nbasis])
        self._S = self.S.ctypes.data_as(POINTER(c_double))

        self.T = np.zeros([self.nbasis, self.nbasis])
        self._T = self.T.ctypes.data_as(POINTER(c_double))

        self.V = np.zeros([self.nbasis, self.nbasis])
        self._V = self.V.ctypes.data_as(POINTER(c_double))

        self.H = np.zeros([self.nbasis, self.nbasis])
        self._H = self.H.ctypes.data_as(POINTER(c_double))

        self.C = np.zeros([self.nbasis, self.nbasis])
        self._C = self.C.ctypes.data_as(POINTER(c_double))

        self.D = np.zeros([self.nbasis, self.nbasis])
        self._D = self.D.ctypes.data_as(POINTER(c_double))

        self.L = np.zeros([self.nbasis, self.nbasis])
        self._L = self.L.ctypes.data_as(POINTER(c_double))

        self.G = np.zeros([self.nbasis, self.nbasis])
        self._G = self.G.ctypes.data_as(POINTER(c_double))

        self.J = np.zeros([self.nbasis, self.nbasis])
        self._J = self.J.ctypes.data_as(POINTER(c_double))

        self.F = np.zeros([self.nbasis, self.nbasis])
        self._F = self.F.ctypes.data_as(POINTER(c_double))

        self.O = np.zeros(self.nbasis)
        self._O = self.O.ctypes.data_as(POINTER(c_double))

        if self.symbol != 'e-beta':
            self._pce = self._compute_pce()

        self.C[:] = psi.C
        self.O[:] = psi.O
        self.D[:] = psi.D

        # Calculate 1 body
        with napmo.runtime.timeblock('Numerical 1 body'):
            self.compute_overlap()
            self.compute_kinetic()
            self.compute_nuclear()
            self.compute_hcore()

        # Calculate initial density
        self.compute_guess()

        # Calculate initial V_el(J)
        self.KO = np.zeros(
            [int(self.nbasis * (self.nbasis + 1) / 2), self._grid.size])

        self.JO = np.zeros(self._grid.size)

        if self.species.get('size') > 1:
            self.KO = self._compute_2body_exchange()

        # # First residual
        # self.R = self._compute_residual()

        # # Build \Delta e (eq. 11)
        # self.DO = self._compute_energy_correction()

        # # Calculate \Delta \phi ^ {(n)}(eq. 13) and Sum \Delta \phi(eq. 12)
        # self._dpsi = self._compute_delta_psi()

        # # Debug information(check eq. 18)
        # print(np.array([np.sqrt(self._grid.integrate(dp * dp))
        #                 for dp in self._dpsi]).sum())

        # Compute \Psi (eq. 14) through conventional SCF for \psi = aC + bDC

        # Calculate \rho = \Psi^2
        # Calculate V^{next}

    def _compute_nuclear_operator(self):
        """
        Computes the action of the Nuclear attraction repulsion
        on a trial function

        :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|} \phi_i`
        """
        op = napmo.compute_nuclear(self._grid, self._pc)

        VO = np.array([phi * op for phi in self.psi])
        VO *= self.species.get('charge')

        # Debug information
        # test_integration(self._grid, self.nbasis,
        #                  self.psi, VO, self.V)

        return VO

    def _compute_kinetic_operator(self):
        """
        Computes the action of the Kinetic operator on a trial function

        :math:`-\dfrac{1}{2} \nabla^{2} \phi_{i}`
        """

        TO = np.array([napmo.compute_kinetic(self._grid, phi, self.lmax)
                       for phi in self.psi])

        # Debug information
        # test_integration(self._grid, self.nbasis,
        #                  self.psi, TO, self.T)

        return TO

    def _build_hcore_operator(self):
        """
        Builds Hamiltonian operator from a previous V_nuc potential and
        Kinetic energy calculation
        """
        HO = np.array([TO + VO for TO, VO in zip(self.TO, self.VO)])

        # Debug information
        # test_integration(self._grid, self.nbasis,
        #                  self.psi, H, self.H)

        return HO

    def _compute_dens_i_from_dm(self):
        """
        Computes the density in the grid for each orbital from density matrix
        """
        with napmo.runtime.timeblock('Numerical Density'):
            DG = np.array([phi * self.D.dot(phi) for phi in self.psi.T]).T

        # Debug  information (Suppose to be the # of e-)
        # print('Number of -e: ', self._grid.integrate(DG.sum(axis=0)))

        return DG

    def _compute_dens_from_dm(self):
        """
        Computes the density in the grid from a density matrix
        """
        DG = np.zeros(self._grid.size)

        # 0.0 -> epsilon
        napmo.cext.nwavefunction_compute_density_from_dm(
            self.species.get('basis')._this,
            self._grid._this,
            self.D, DG, 0.0, np.abs(self.D).max(axis=0))

        # Debug  information (Suppose to be the # of e-)
        # print('Number of -e: ', self._grid.integrate(DG))

        return DG

    def _compute_2body_exchange(self):
        """
        Computes the exchange potential solving Poisson's equation
        following the procedure of Shiozaki, T., & Hirata, S. (2007).
        Grid-based numerical Hartree-Fock solutions of polyatomic molecules.
        Physical Review A - Atomic, Molecular, and Optical Physics, 76(4), 040503.
        http://doi.org/10.1103/PhysRevA.76.040503
        """
        with napmo.runtime.timeblock('Numerical exchange'):
            KO = np.array([napmo.compute_coulomb(
                self._grid, self.psi[s, :] * self.psi[v, :], self.lmax)
                for s in range(self.nbasis)
                for v in range(self.nbasis)
                if s >= v])

        return KO

    def _compute_2body_coulomb(self):
        """
        Computes coulomb potential solving Poisson's equation
        following Becke procedure using the current density
        """
        with napmo.runtime.timeblock('Numerical coulomb'):
            JO = napmo.compute_coulomb(
                self._grid, self.DG.sum(axis=0), self.lmax)

        # Debug information
        # print("Coulomb energy: ", 0.5 *
        #       self._grid.integrate(J * self.D.sum(axis=0)))

        return JO

    def _compute_residual(self):
        """
        Build R (eq. 10) :math:`R = (T + V - e) \phi`
        """
        O = np.zeros([self.nbasis, self.nbasis])
        np.fill_diagonal(O, self.O)

        aux = self.F.dot(self.C) - self.S.dot(self.C.dot(O))
        R = np.array([aux.dot(phi) for phi in self.psi.T]).T

        # print(aux.sum())

        return R

    def _compute_energy_correction(self):
        """
        Computes energy correction. eq. 11
        """
        DO = np.array([self._grid.integrate(r * phi)
                       for r, phi in zip(self.R, self.psi)])

        # print(DO.sum())

        return DO

    def get_operator_matrix(self, Op):
        """
        Calculates the matrix representation of the operator ``Op``
        based on the grid representation of such operator.

        Each element on the matrix corresponds to:

        :math:`\int \psi_i Op \psi_j dr^3`

        Args:
            O (ndarray): :math:`Op \psi` calculated on the grid.
        """
        O = np.array([self._grid.integrate(self.psi[i, :] * Op[j])
                      for i in range(self.nbasis)
                      for j in range(self.nbasis)])

        O = O.reshape([self.nbasis, self.nbasis])

        return O

    def _compute_delta_psi(self):
        """
        Computes \Delta \psi. eq. 13
        """
        dpsi = np.array([napmo.compute_dpsi(
            self._grid, self.lmax, phi, doi, oi, ri, self.JO + vnuc)
            for phi, doi, oi, ri, vnuc in zip(
                self.psi, self.DO, self.O, self.R, self.VO)])

        return dpsi

    def compute_overlap(self):
        """
        Computes the overlap matrix
        """
        self.S[:] = self.get_operator_matrix(self.psi)
        # print("\n Overlap Matrix:")
        # print(self.S)

    def compute_kinetic(self):
        """
        Computes the Kinetic matrix
        """
        self.TO = self._compute_kinetic_operator()
        self.T[:] = self.get_operator_matrix(self.TO)
        # print("\n Kinetic Matrix:")
        # print(self.T)

    def compute_nuclear(self):
        """
        Computes the nuclear-electron potential matrix
        """
        self.VO = self._compute_nuclear_operator()
        self.V[:] = self.get_operator_matrix(self.VO)
        # print("\n Attraction Matrix")
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix
        """

        if self.species.get('size') > 1:

            self.DG = self._compute_dens_i_from_dm()
            self.JO = self._compute_2body_coulomb()

            with napmo.runtime.timeblock('Numerical G matrix'):
                napmo.cext.nwavefunction_compute_2body_matrix(
                    byref(self), self._grid._this, self.psi, self.JO, self.KO)

            self.G *= self.species.get('charge')**2

        # print("\n Coupling Matrix:")
        # print(self.J)

    def compute_hcore(self):
        """
        Builds the Hcore matrix
        """
        self.HO = self._build_hcore_operator()
        self.H[:] = self.T + self.V
        # print("\n hcore Matrix")
        # print(self.H)

    def compute_guess(self):
        """
        Computes the density guess using analytical density.
        """
        self.DG = self._compute_dens_i_from_dm()

    def build_fock(self):
        """
        Builds the Fock matrix
        """

        self.F[:] = self.H + self.G + self.J
        # print("\n Fock Matrix:")
        # print(self.F)

    def _compute_pce(self):
        """
        Calculates the total point charges energy for the system
        """
        output = [pa.get('charge') * pb.get('charge') /
                  np.sqrt(((pa.get('origin') - pb.get('origin'))**2).sum())
                  for i, pa in enumerate(self.species.get('particles'))
                  for j, pb in enumerate(self.species.get('particles'))
                  if j > i and not pa.is_quantum and not pb.is_quantum
                  ]
        return sum(output)

    @property
    def nbasis(self):
        """
        Length of the basis
        """
        return self._nbasis

    @property
    def lmax(self):
        """
        Maximum ``l`` for the spherical expansion
        """
        return self._lmax

    @property
    def pce(self):
        """
        The potential charge energy
        """
        return self._pce

    @property
    def symbol(self):
        """
        Symbol of the object owner's
        """
        return self._symbol

    @property
    def sid(self):
        """
        Id of the object owner's
        """
        return self._sid
