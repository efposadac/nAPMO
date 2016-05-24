# file: libint.py
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

from napmo.system.cext import napmo_library as nl
from napmo.system.basis_set import BasisSet_C


class Libint(object):

    """
    Libint Library interface
    """

    def __init__(self, molecule):
        super(Libint, self).__init__()

        self._this = nl.LibintInterface_new()

        for atom in molecule.get('atoms'):
            print(molecule.get_basis(atom.get('symbol')))
            nl.LibintInterface_add_particle(
                self._this,
                atom.get('atomic_number'),
                atom.get('origin'),
                byref(molecule.get_basis_as_cstruct(atom.get('symbol'))))

    def __del__(self):
        nl.LibintInterface_del(self._this)

    def get_overlap_matrix(self):
        nbasis = nl.LibintInterface_get_nbasis(self._this)
        output = np.empty([nbasis, nbasis])
        nl.LibintInterface_compute_1body_ints(self._this, 1, output)
        return output

    def get_kinetic_matrix(self):
        nbasis = nl.LibintInterface_get_nbasis(self._this)
        output = np.empty([nbasis, nbasis])
        nl.LibintInterface_compute_1body_ints(self._this, 2, output)
        return output

    def get_nuclear_matrix(self):
        nbasis = nl.LibintInterface_get_nbasis(self._this)
        output = np.empty([nbasis, nbasis])
        nl.LibintInterface_compute_1body_ints(self._this, 3, output)
        return output


array_1d_double = npct.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array_2d_double = npct.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

nl.LibintInterface_new.restype = c_void_p

nl.LibintInterface_del.restype = None

nl.LibintInterface_add_particle.restype = None
nl.LibintInterface_add_particle.argtypes = [
    c_void_p, c_int, array_1d_double, POINTER(BasisSet_C)
]

nl.LibintInterface_compute_1body_ints.restype = None
nl.LibintInterface_compute_1body_ints.argtypes = [
    c_void_p, c_int, array_2d_double]

nl.LibintInterface_get_nbasis.restype = c_int
nl.LibintInterface_get_nbasis.argtypes = [c_void_p]
