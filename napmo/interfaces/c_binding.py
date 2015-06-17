# file: c_binding.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from ctypes import *
from copy import deepcopy
import os
import sysconfig

from napmo.utilities.constants import *


class C_BasisSet(Structure):
    """
    C interface to the BasisSet class.
    """
    _fields_ = [
        ("n_cont", c_int),
        ("n_prim_cont", POINTER(c_int)),
        ("l_index", POINTER(c_int)),
        ("exponent", POINTER(c_double)),
        ("coefficient", POINTER(c_double)),  # normalized
        ("normalization", POINTER(c_double)),
        ("origin", POINTER(c_double))
    ]

    def __init__(self, basis):
        super(C_BasisSet, self).__init__()

        self.n_cont = basis.get('length')

        size = (c_double * self.n_cont)()
        self.normalization = cast(size, POINTER(c_double))

        size = (c_double * self.n_cont * 3)()
        self.origin = cast(size, POINTER(c_double))

        size = (c_int * self.n_cont)()
        self.n_prim_cont = cast(size, POINTER(c_int))

        size = (c_int * self.n_cont * 3)()
        self.l_index = cast(size, POINTER(c_int))

        size = (c_double * basis.get('t_length'))()
        self.exponent = cast(size, POINTER(c_double))

        size = (c_double * basis.get('t_length'))()
        self.coefficient = cast(size, POINTER(c_double))

        counter = 0
        for i in range(self.n_cont):
            F = basis.get('function')[i]
            self.normalization[i] = F.get('normalization')

            for j in range(3):
                self.l_index[i*3+j] = F.get('l')[j]
                self.origin[i*3+j] = F.get('origin')[j]

            self.n_prim_cont[i] = F.get('length')
            for j in range(self.n_prim_cont[i]):
                P = F.get('primitive')[j]
                self.exponent[counter] = P.get('exponent')
                self.coefficient[counter] = P.get('coefficient') * P.get('normalization')
                counter += 1


class CBinding(Structure):
    """
    Interfaces to C napmo library.
    """
    _fields_ = [
        ("n_particles", c_int),
        ("particle_number", POINTER(c_int)),
        ("particle_radii", POINTER(c_double)),
        ("particle_origin", POINTER(c_double)),
        ("work_space", POINTER(c_double)),
        ("basis_set", C_BasisSet)
    ]

    def __init__(self, particle_stack):
        super(CBinding, self).__init__()

        # System information
        self.n_particles = len(particle_stack)

        # Particle information
        size = (c_double * self.n_particles)()
        self.particle_radii = cast(size, POINTER(c_double))

        size = (c_int * self.n_particles)()
        self.particle_number = cast(size, POINTER(c_int))

        size = (c_double * self.n_particles * 3)()
        self.particle_origin = cast(size, POINTER(c_double))

        size = (c_double * self.n_particles)()
        self.work_space = cast(size, POINTER(c_double))

        basis = deepcopy(particle_stack[0].get('basis'))

        for i in range(self.n_particles):
            atom = particle_stack[i].get
            self.particle_radii[i] = atom('atomic_radii') * ANGSTROM_TO_BOHR
            self.particle_number[i] = atom('atomic_number')

            for j in range(3):
                self.particle_origin[i*3+j] = atom('origin')[j]

            if i > 0:
                basis += atom('basis')

        # Basis-set information
        self.basis_set = C_BasisSet(basis)

#################################################
# Interface to napmo_library
#################################################

try:
    lp = os.path.join(os.path.dirname(__file__), '../../libnapmo'+sysconfig.get_config_var('SO'))
    napmo_library = CDLL(lp)
except OSError:
    lp = os.path.join(os.path.dirname(__file__), '../../src/libnapmo.so')
    napmo_library = CDLL(lp)
