# file: basis_set.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import json
from ctypes import *

from napmo.system.contracted_slater import ContractedSlater
from napmo.system.contracted_gaussian import ContractedGaussian


class BasisSet(dict):

    """
    Basis-set interface. (dict)

    This class allows the management of STO and GTO basis-sets.
    """

    def __init__(self, name='user'):
        super(BasisSet, self).__init__()
        self['name'] = name
        self['function'] = []
        self['kind'] = None
        self['length'] = 0
        self['t_length'] = 0

    def load_gaussian(self, particle, data,
                      origin=np.zeros(3, dtype=np.float64)):
        """
        Load a Gaussian Type Orbital (GTO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self['kind'] = 'GTO'
        data = json.loads(data)[particle]

        lvalue = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

        func = [ContractedGaussian(
            np.array(cont['prim'], dtype=np.float64),
            np.array(cont['cont'], dtype=np.float64),
            np.array(origin, dtype=np.float64),
            np.array([i, k-j, j], dtype=np.int32))
            for cont in data
            for k, i in enumerate(reversed(
                range(lvalue.get(cont['angular']) + 1)))
            for j in range(k + 1)]

        self['function'] += func
        self['length'] = len(self.get('function'))
        self['t_length'] = 0

        self['t_length'] = sum([function.get('length')
                                for function in self.get('function')])

    def load_slater(self, particle, data,
                    origin=np.zeros(3, dtype=np.float64)):
        """
        Load a Slater Type Orbital (STO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self["kind"] = 'STO'
        data = json.loads(data)[particle]

        lvalue = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}

        func = [ContractedSlater(
            np.array(cont["prim"], dtype=np.float64),
            np.array(cont["cont"], dtype=np.float64),
            np.array(origin, dtype=np.float64),
            np.array(cont["n"], dtype=np.int32),
            lvalue.get(cont['angular']),
            m)
            for cont in data
            for m in range(-lvalue.get(cont['angular']),
                           lvalue.get(cont['angular']) + 1)]

        self['function'] += func
        self['length'] = len(self.get('function'))
        self['t_length'] = 0

        self['t_length'] = sum([function.get('length')
                                for function in self.get('function')])

    def compute(self, coord=np.zeros([1, 3], dtype=np.float64)):
        """
        Compute all the basis-set functions at given ``coord``.

        Args:
            coord (ndarray): coordinates in which the basis set will be
                evaluated. Array shape should be (n, 3)

        """
        return [func.compute(coord) for func in self.get('function')]

    def __repr__(self):
        """
        Prints extended information of the basis-set object.
        """
        out = ('================\n' +
               'Basis set info  \n' +
               '================\n' +
               'Name: '+self.get('name') + '\n'
               'Kind: '+self.get('kind') + '\n'
               'Length: '+str(self.get('length')) + '\n'
               '****************\n' +
               'Functions info: \n' +
               '****************')

        out += "\n".join([str(p) for p in self.get('function')])

        return out


class BasisSet_C(Structure):

    """
    C interface to the BasisSet class.

    Args:
        basis (BasisSet) : Basis set of the system.
    """
    _fields_ = [
        ("n_cont", c_int),
        ("n_prim_cont", POINTER(c_int)),
        ("prim_index", POINTER(c_int)),
        ("l_index", POINTER(c_int)),
        ("origin", POINTER(c_double)),
        ("normalization", POINTER(c_double)),
        ("exponent", POINTER(c_double)),
        ("coefficient", POINTER(c_double))  # normalized
    ]

    def __init__(self, basis):
        super(BasisSet_C, self).__init__()

        self.n_cont = basis.get('length')

        size = (c_int * self.n_cont)()
        self.n_prim_cont = cast(size, POINTER(c_int))

        size = (c_int * self.n_cont)()
        self.prim_index = cast(size, POINTER(c_int))

        size = (c_int * self.n_cont * 3)()
        self.l_index = cast(size, POINTER(c_int))

        size = (c_double * self.n_cont * 3)()
        self.origin = cast(size, POINTER(c_double))

        size = (c_double * self.n_cont)()
        self.normalization = cast(size, POINTER(c_double))

        size = (c_double * basis.get('t_length'))()
        self.exponent = cast(size, POINTER(c_double))

        size = (c_double * basis.get('t_length'))()
        self.coefficient = cast(size, POINTER(c_double))

        counter = 0
        for i in range(self.n_cont):
            F = basis.get('function')[i]
            self.normalization[i] = F.get('normalization')

            for j in range(3):
                self.l_index[i * 3 + j] = F.get('l')[j]
                self.origin[i * 3 + j] = F.get('origin')[j]

            self.n_prim_cont[i] = F.get('length')
            self.prim_index[i] = counter

            for j in range(self.n_prim_cont[i]):
                P = F.get('primitive')[j]
                self.exponent[counter] = P.exponent
                self.coefficient[counter] = P.coefficient * P.normalization
                counter += 1
