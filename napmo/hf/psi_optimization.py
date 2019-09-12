from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
from scipy.optimize import minimize

import napmo


def INDEX(i, j):
    if i > j:
        return int((i * (i + 1) / 2) + j)
    else:
        return int((j * (j + 1) / 2) + i)


class PSIO(napmo.PSIB):
    """
    Class for wavefunction optimization. Solves :math:`\Psi = a \Psi + \Delta \Psi` using a N2-Dim
    linear variational optimization
    """

    def __init__(self, psin):
        super(PSIO, self).__init__(psin.species,
                                   ndim=psin.ndim * 2)

        self._debug = psin._debug
        self._diis = napmo.cext.LibintInterface_diis_new(2)
        self._pc = psin._pc
        self.species = psin.species
        self._grid = psin._grid
        self._lmax = psin._lmax
        self.Vnuc = psin.Vnuc
        self._mass_inv = psin._mass_inv

        self._exchange = False
        self.iterations = 0

        self.Kgrid = np.zeros([self.ndim, self._grid.size])
        self.Jgrid = np.zeros(self._grid.size)
        self.XCgrid = np.zeros(self._grid.size)

    def optimize(self, wf, scf, other_wf=None):
        """
        Calculates \Delta \phi as Eq. 14 Becke's paper.
        """
        self._exchange = False
        self._rmsd = 1.0
        self.psi = np.vstack([wf.psi, wf.delta_psi])

        self.normalize()
        self.compute_1body()

        if self.iterations == 0:
            self.compute_guess()

        scf.single(self, pprint=self._debug, diis=False, other_psi=other_wf)
        self.iterations += 1

        return self._compute_psi_from_cm(self.C, self.psi), self.O[:self._occupation]

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
        self.V *= self.species.get('charge')

        # print("\n Attraction Matrix")
        # print(self.V)

    def compute_2body(self, direct=False):
        """
        Computes the two-body matrix.
        """

        self._compute_density()

        if self.species.get('size') > 1:
            self._compute_2body_coulomb()
            if self._exchangefactor != 0.0:
                self._compute_2body_exchange()
            #FELIX: add if bla bla bla bla
            # self._compute_2body_exchange()
            with napmo.runtime.timeblock('Numerical G matrix'):
                napmo.cext.nwavefunction_compute_2body_matrix_atm(
                byref(self), self._grid._this, self.psi, self.Jgrid, self.Kgrid)

            self.G *= self.species.get('charge')

        # print("\n G Matrix:")
        # print(self.G)

    def compute_coupling(self, other_psi, direct=False):
        """
        Computes the two-body coupling matrix

        Args:
            other_psi (WaveFunction) : WaveFunction object for the other species.
            direct (bool) : Whether to calculate eris on-the-fly or not
        """
        aux = np.zeros([self._ndim, self._ndim])
        self.J[:] = 0.0
        for psi in other_psi:
            if self.sid != psi.sid:

                psi.Cgrid = napmo.compute_coulomb(
                    psi._grid, psi.Dgrid.sum(axis=0), psi.lmax)

                psi.Cgrid *= self.species.get('charge') * \
                    psi.species.get('charge')

                napmo.cext.nwavefunction_compute_coupling(
                    byref(self), self._grid._this, self.psi, psi.Cgrid, aux)

                self.J += aux

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

        self.F[:] = self.H + self.G + self.J + self.XC

        # print("\n Fock Matrix:")
        # print(self.F)

    def _compute_density(self):
        """
        Compute the density on the grid
        """

        self.Dgrid = self._compute_density_from_dm(
            self.D, self.psi)

        # print("\n Density Matrix: ", self.symbol)
        # print(self.D)

        # Debug information (Suppose to be the # of e-)
        # print('Number of -e: ', self._grid.integrate(self.Dgrid.sum(axis=0)))

    def _compute_density_from_dm(self, dm, psi):
        """
        Computes the density in the grid for each orbital from density matrix. Includes virtual orbitals.
        """
        with napmo.runtime.timeblock('Numerical Density'):
            dens = np.array([phi * dm.dot(phi)
                             for phi in psi.T]).T

        return dens

    def _compute_psi_from_cm(self, cm, psi):

        res = np.zeros([self._occupation, self._grid.size])

        with napmo.runtime.timeblock('Numerical Density'):
            for i in range(self._occupation):
                for j in range(cm.shape[0]):
                    res[i] += cm[j, i] * psi[j]

        # Debug information(Suppose to be the  # of e-)
        # print('Number of -e: ', self._grid.integrate(res.sum(axis=0)
        #                                              * res.sum(axis=0)) * self._eta)

        return res

    def _compute_nuclear_operator(self):
        """
        Computes the action of the Nuclear attraction repulsion
        on a trial function

        :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|} \phi_i`

        """

        self.Vgrid = np.array([phi * self.Vnuc for phi in self.psi])

    def _compute_kinetic_operator(self):
        """
        Computes the action of the Kinetic operator on a trial function

        :math:`-\dfrac{1}{2} \nabla^{2} \phi_{i}`
        """

        self.Tgrid = np.array([napmo.compute_kinetic(self._grid, phi, self.lmax)
                               for phi in self.psi])

        self.Tgrid *= self._mass_inv

    def _compute_2body_coulomb(self):
        """
        Computes coulomb potential solving Poisson's equation
        following Becke procedure using the current density
        """
        if self.species.get('size') > 1:

            with napmo.runtime.timeblock('Numerical coulomb'):

                self.Jgrid[:] = napmo.compute_coulomb(
                    self._grid, self.Dgrid.sum(axis=0), self.lmax)

                self.Jgrid *= self.species.get('charge')

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

                    aux = np.array([napmo.compute_coulomb(
                        self._grid, self.psi[s, :] * self.psi[v, :], self.lmax)
                        for s in range(self.ndim)
                        for v in range(self.ndim)
                        if s >= v])

                    self.Kgrid[:] = np.array([np.array([self.psi[k, :] * aux[int(INDEX(l, j))] * self.D[k, l]
                                                        for k in range(self.ndim)
                                                        for l in range(self.ndim)]).sum(axis=0)
                                              for j in range(self.ndim)])

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
        # print ("Trololo se cae aqui", self.ndim, "size O", O.size)

        M = np.array([self._grid.integrate(self.psi[i] * O[j])
                      for i in range(self.ndim)
                      for j in range(self.ndim)])

        M = M.reshape([self.ndim, self.ndim])

        # M = np.zeros([self.ndim, self.ndim])
        # for i in range(self.ndim):
        #     for j in range(self.ndim):
        #         M[i,j] = self._grid.integrate(self.psi[i] * O[j])
        
        return M

    @property
    def lmax(self):
        return self._lmax
