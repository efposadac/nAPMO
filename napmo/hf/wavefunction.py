# file: wavefunction.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import numpy.ctypeslib as npct

import napmo


class WaveFunction(Structure):

    """
    """
    _fields_ = [
        ("_S", POINTER(c_double)),
        ("_T", POINTER(c_double)),
        ("_V", POINTER(c_double)),
        ("_H", POINTER(c_double)),
        ("_C", POINTER(c_double)),
        ("_D", POINTER(c_double)),
        ("_L", POINTER(c_double)),
        ("_G", POINTER(c_double)),
        ("_J", POINTER(c_double)),
        ("_F", POINTER(c_double)),
        ("_nbasis", c_int),
        ("_occupation", c_int),
        ("_eta", c_double),
        ("_kappa", c_double),
        ("_energy", c_double),
        ("_rmsd", c_double)
    ]

    def __init__(self, species, point_charges):
        super(WaveFunction, self).__init__()
        self.species = species

        # Initialize Libint object to calculate integrals
        self._libint = napmo.cext.LibintInterface_new(species.get('id'))

        for particle in species.get('particles'):
            napmo.cext.LibintInterface_add_basis(
                self._libint,
                byref(napmo.BasisSet_C(particle.get('basis'))))

        for point in point_charges:
            napmo.cext.LibintInterface_add_pointcharges(
                self._libint,
                point.get('atomic_number', point.get('charge')),
                point.get('origin'))

        # Initialize defaults
        self._nbasis = napmo.cext.LibintInterface_get_nbasis(self._libint)
        self._occupation = species.get('occupation')
        self._eta = species.get('eta')
        self._kappa = species.get('kappa')
        self._ints = None
        self._symbol = species.get('symbol')
        self._sid = species.get('id')
        self._energy = 0.0
        self._rmsd = 1.0
        self._pce = 0.0

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

        if self.symbol != 'e-beta':
            self._pce = self._compute_pce()

    def compute_overlap(self):
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 1, self.S)
        # print("\n Overlap Matrix:")
        # print(self.S)

    def compute_kinetic(self):
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 2, self.T)
        self.T /= self.species.get('mass')
        # print("\n Kinetic Matrix:")
        # print(self.T)

    def compute_nuclear(self):
        napmo.cext.LibintInterface_compute_1body_ints(self._libint, 3, self.V)
        self.V *= -self.species.get('charge')
        # print("\n Attraction Matrix")
        # print(self.V)

    def compute_2body(self, direct=False):
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

        # print("\n G Matrix:")
        # print(self.G)

    def compute_couling(self, other_psi, direct=False):
        aux = np.zeros([self.nbasis, self.nbasis])
        self.J[:] = 0.0
        for psi in other_psi:
            if self.sid != psi.sid:
                napmo.cext.LibintInterface_compute_coupling_direct(
                    self._libint, psi._libint, psi.D, aux)
                aux *= self.species.get('charge') * psi.species.get('charge')
                self.J += aux

        self.J += np.tril(self.J, -1).T

        # print("\n Coupling Matrix:")
        # print(self.J)

    def compute_hamiltonian(self):
        self.H[:] = self.T + self.V
        # print("\n Hamiltonian Matrix")
        # print(self.H)

    def compute_guess(self):
        napmo.cext.wavefunction_guess_hcore(byref(self))
        self.D *= self._eta
        # print("\n Coefficients Matrix")
        # print(self.C)

        # print("\n Density Guess")
        # print(self.D)

    def build_fock(self):
        self.F[:] = self.H + self.G + self.J
        # print("\n Fock Matrix:")
        # print(self.F)

    def _compute_pce(self):

        output = [pa.get('charge') * pb.get('charge') /
                  np.sqrt(((pa.get('origin') - pb.get('origin'))**2).sum())
                  for i, pa in enumerate(self.species.get('particles'))
                  for j, pb in enumerate(self.species.get('particles'))
                  if j > i and not pa.is_quantum and not pb.is_quantum
                  ]

        return sum(output)

    @property
    def nbasis(self):
        return self._nbasis

    @property
    def pce(self):
        return self._pce

    @property
    def symbol(self):
        return self._symbol

    @property
    def sid(self):
        return self._sid


array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo.cext.LibintInterface_new.restype = c_void_p
napmo.cext.LibintInterface_new.argtypes = [c_int]

napmo.cext.LibintInterface_del.restype = None

napmo.cext.LibintInterface_add_pointcharges.restype = None
napmo.cext.LibintInterface_add_pointcharges.argtypes = [
    c_void_p, c_int, array_1d_double]

napmo.cext.LibintInterface_add_basis.restype = None
napmo.cext.LibintInterface_add_basis.argtypes = [
    c_void_p, POINTER(napmo.BasisSet_C)]

napmo.cext.LibintInterface_compute_1body_ints.restype = None
napmo.cext.LibintInterface_compute_1body_ints.argtypes = [
    c_void_p, c_int, array_2d_double]

napmo.cext.LibintInterface_get_nbasis.restype = c_int
napmo.cext.LibintInterface_get_nbasis.argtypes = [c_void_p]

napmo.cext.LibintInterface_init_2body_ints.restype = None
napmo.cext.LibintInterface_init_2body_ints.argtypes = [c_void_p]

napmo.cext.LibintInterface_compute_2body_ints.restype = c_void_p
napmo.cext.LibintInterface_compute_2body_ints.argtypes = [
    c_void_p, array_2d_double]

napmo.cext.LibintInterface_compute_2body_direct.restype = None
napmo.cext.LibintInterface_compute_2body_direct.argtypes = [
    c_void_p, array_2d_double, array_2d_double]

napmo.cext.LibintInterface_compute_coupling_direct.restype = None
napmo.cext.LibintInterface_compute_coupling_direct.argtypes = [
    c_void_p, c_void_p, array_2d_double, array_2d_double]

napmo.cext.wavefunction_guess_hcore.restype = None
napmo.cext.wavefunction_guess_hcore.argtypes = [POINTER(WaveFunction)]

napmo.cext.wavefunction_compute_2body_matrix.restype = None
napmo.cext.wavefunction_compute_2body_matrix.argtypes = [
    POINTER(WaveFunction), c_void_p]
