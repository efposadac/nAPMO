# file: psi_base.py
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


class PSIB(Structure):
    """
    Defines the Fock operator for a Hartree-Fock Calculation.
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
        ("_ndim", c_int),
        ("_occupation", c_int),
        ("_eta", c_double),
        ("_kappa", c_double),
        ("_energy", c_double),
        ("_rmsd", c_double)  # Root medium square deviation, for D matrix
    ]

    def __init__(self, species, ndim=None):
        super(PSIB, self).__init__()
        self.species = species

        self._nbasis = species.get('basis').get('nbasis')
        self._occupation = species.get('occupation')
        self._eta = species.get('eta')
        self._kappa = species.get('kappa')
        self._ints = None
        self._symbol = species.get('symbol')
        self._sid = species.get('id')
        self._energy = 0.0
        self._rmsd = 1.0
        self._pce = 0.0
        self._energy = 0.0
        self._pce = 0.0

        if ndim is None:
            self._ndim = self.nbasis
        else:
            self._ndim = ndim

        self.S = np.zeros([self._ndim, self._ndim])
        self._S = self.S.ctypes.data_as(POINTER(c_double))

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
