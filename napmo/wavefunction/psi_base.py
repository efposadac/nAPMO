# file: psi_base.py
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

# import matplotlib.pyplot as plt


class PSIB(Structure):
    """
    Defines the Fock operator for a Hartree-Fock Calculation.
    """
    _fields_ = [
        ("_S", POINTER(c_double)),   # Overlap
        ("_X", POINTER(c_double)),   # Transformation
        ("_T", POINTER(c_double)),   # Kinetic
        ("_V", POINTER(c_double)),   # Nuclear
        ("_H", POINTER(c_double)),   # Hcore
        ("_C", POINTER(c_double)),   # Coefficients
        ("_D", POINTER(c_double)),   # Density
        ("_L", POINTER(c_double)),   # Last Density
        ("_G", POINTER(c_double)),   # 2 Body
        ("_J", POINTER(c_double)),   # Coupling
        ("_XC", POINTER(c_double)),  # Exchange Correlation Matrix
        ("_F", POINTER(c_double)),   # Fock
        ("_O", POINTER(c_double)),   # Orbitals
        ("_nbasis", c_int),
        ("_ndim", c_int),
        ("_occupation", c_int),
        ("_eta", c_double),
        ("_kappa", c_double),
        ("_x_factor", c_double),      # Fraction of exchange
        ("_xc_energy", c_double),     # Exchange Correlation Energy
        ("_energy", c_double),
        ("_rmsd", c_double)           # Root medium square deviation, for D matrix
    ]

    def __init__(self, species, ndim=None, options={}):
        super(PSIB, self).__init__()
        self.species = species

        self.options = {}
        self.options.update(options)

        self._nbasis = species.get('basis').get('nbasis')
        self._occupation = species.get('occupation')
        self._eta = species.get('eta')
        self._kappa = species.get('kappa')
        self._functional = None
        self._ints = None
        self._symbol = species.get('symbol')
        self._sid = species.get('id')
        self._rmsd = 1.0
        self._pce = 0.0
        self._energy = 0.0
        self._xc_energy = 0.0
        self._xc_vrho = None
        self._tf = False

        self._x_factor = self._kappa / self._eta

        if 'tf' in self.options:
            self._tf = True
            print("Using Translation-Free Correction!!!")

        self._method = self.options.get('method')

        if self._method == 'dft':
            if self.species.get('is_electron') and self.symbol != 'e-beta':
                self._functional = napmo.Functional(self._symbol, self.options)
                # TODO: Fix this
                # self._x_factor = self._functional.x_factor
            self._x_factor = 0.0
            # self._x_factor = self._kappa/self._eta

        if ndim is None:
            self._ndim = self.nbasis
        else:
            self._ndim = ndim

        assert self._ndim > 0

        self.S = np.zeros([self._ndim, self._ndim])
        self._S = self.S.ctypes.data_as(POINTER(c_double))

        self.X = np.zeros([self._ndim, self._ndim])
        self._X = self.X.ctypes.data_as(POINTER(c_double))

        self.T = np.zeros([self._ndim, self._ndim])
        self._T = self.T.ctypes.data_as(POINTER(c_double))

        self.V = np.zeros([self._ndim, self._ndim])
        self._V = self.V.ctypes.data_as(POINTER(c_double))

        self.H = np.zeros([self._ndim, self._ndim])
        self._H = self.H.ctypes.data_as(POINTER(c_double))

        self.C = np.zeros([self._ndim, self._ndim])
        self._C = self.C.ctypes.data_as(POINTER(c_double))

        self.D = np.zeros([self._ndim, self._ndim])
        self._D = self.D.ctypes.data_as(POINTER(c_double))

        self.L = np.zeros([self._ndim, self._ndim])
        self._L = self.L.ctypes.data_as(POINTER(c_double))

        self.G = np.zeros([self._ndim, self._ndim])
        self._G = self.G.ctypes.data_as(POINTER(c_double))

        self.J = np.zeros([self._ndim, self._ndim])
        self._J = self.J.ctypes.data_as(POINTER(c_double))

        self.XC = np.zeros([self._ndim, self._ndim])
        self._XC = self.XC.ctypes.data_as(POINTER(c_double))

        self.F = np.zeros([self._ndim, self._ndim])
        self._F = self.F.ctypes.data_as(POINTER(c_double))

        self.O = np.zeros(self._ndim)
        self._O = self.O.ctypes.data_as(POINTER(c_double))

        if self.symbol != 'e-beta':
            self._pce = self._compute_pce()

    def _compute_pce(self):
        """
        Calculates the total point charges energy for the system
        """
        output = [
            pa.get('charge') * pb.get('charge') /
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
    def ndim(self):
        return self._ndim

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

    @property
    def energy(self):
        return self._energy

    @property
    def xc_energy(self):
        return self._xc_energy

    @property
    def occupation(self):
        return self._occupation

    # ### TEST: DO NOT use this functions for real things ###

    def _compute_density_from_dm(self, dm, psi):
        """
        Computes the density in the grid for each orbital from density matrix.
        Includes virtual orbitals.
        """
        with napmo.runtime.timeblock('Numerical Density'):
            dens = np.array([phi * dm.dot(phi)
                             for phi in psi.T]).T
        return dens

    def get_data(self, grid, index, F):

        data = np.zeros([index.size, 2])

        for i, j in enumerate(index):
            data[i, 0] = grid.points[j, 2]
            data[i, 1] = F[j]

        return data

    def get_index_center(self):

        index = [np.array([i for i in range(atgrid.size) if np.abs(
            atgrid.points[i, 0]) < 1.0e-10 and np.abs(atgrid.points[i, 1]) < 1.0e-10])
            for atgrid in self._grid.atgrids]

        return index

    def plot_dens(self, psi=None, grid=None, kind='', marker="-", xlim=None, ylim=None):
        print("Plotting Density...", end=' ')

        if grid is not None:
            self._grid = grid

        if psi is None:
            gbasis = self.species.get('basis').compute(self._grid.points).T.copy()
            Dgrid = self._compute_density_from_dm(self.D, gbasis)
        else:
            Dgrid = psi**2 * self._eta

        self.plot_obj(Dgrid.sum(axis=0), label="rho" + self.symbol + kind, marker=marker, xlim=xlim, ylim=ylim)

        print("Done!")

    def plot_obj(self, obj, label, marker="-", xlim=None, ylim=None):
        """
        Computes the density in the grid for each orbital from density matrix.
        Includes virtual orbitals.
        """

        index = self.get_index_center()
        bwa = np.zeros([0, 2])

        for i, atgrid in enumerate(self._grid.atgrids):
            tmp = self.get_data(atgrid, index[i], obj[
                i * atgrid.size:i * atgrid.size + atgrid.size])
            bwa = np.concatenate((bwa, tmp))

        bwa.view("float64, float64").sort(order=["f0"], axis=0)

        np.savetxt(
            label + '.txt', np.stack([bwa[:, 0], bwa[:, 1]]).T, delimiter='\t', fmt=['%.12f', '%.12f'])

        # plt.plot(bwa[:, 0] * napmo.BOHR_TO_ANGSTROM,
        #          bwa[:, 1], marker, label=label)
        # plt.xlabel("z")
        # plt.ylabel("obj")
        # plt.legend()

        # if xlim is not None:
        #     plt.xlim(xlim[0], xlim[1])

        # if ylim is not None:
        #     plt.ylim(ylim[0], ylim[1])
