# file: basis_set.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

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

    def load_gaussian(self, particle, data, origin=np.zeros(3, dtype=np.float64)):
        """
        Load a Gaussian Type Orbital (GTO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self['kind'] = 'GTO'
        data = json.loads(data)[particle]

        lvalue = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

        for k in range(len(data)):
            l = lvalue[data[k]['angular']]
            for i in range(l + 1):
                x = l - i
                for j in range(i + 1):
                    y = i - j
                    z = j
                    self['function'].append(ContractedGaussian(
                        np.array(data[k]['prim'], dtype=np.float64),
                        np.array(data[k]['cont'], dtype=np.float64),
                        np.array(origin, dtype=np.float64),
                        np.array([x, y, z], dtype=np.int32)
                    ))

        self['length'] = len(self.get('function'))
        self['t_length'] = 0

        for function in self.get('function'):
            self['t_length'] += function.get('length')

    def load_slater(self, particle, data, origin=np.zeros(3, dtype=np.float64)):
        """
        Load a Slater Type Orbital (STO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self["kind"] = 'STO'
        data = json.loads(data)[particle]

        lvalue = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}

        for k in range(len(data)):
            l = lvalue[data[k]["angular"]]
            for m in range(-l, l + 1):
                self["function"].append(ContractedSlater(
                    np.array(data[k]["prim"], dtype=np.float64),
                    np.array(data[k]["cont"], dtype=np.float64),
                    np.array(origin, dtype=np.float64),
                    np.array(data[k]["n"], dtype=np.int32),
                    l,
                    m
                ))
        self['length'] = len(self.get('function'))
        self['t_length'] = 0

        for function in self.get('function'):
            self['t_length'] += function.get('length')

    def compute(self, coord=np.zeros([1, 3], dtype=np.float64)):
        """
        Compute all the basis-set functions at given ``coord``.

        Args:
            coord (numpy.ndarray(3)): coordinates in which the basis set will be evaluated.

        """
        output = []
        for i in range(self.get('length')):
            output.append(self.get('function')[i].compute(coord))

        return output

    def show(self):
        """
        Prints extended information of the basis-set object.
        """
        print("================")
        print("Basis set info")
        print("================")
        print("Name: ", self.get('name'))
        print("Kind: ", self.get('kind'))
        print("Length: ", self.get('length'))
        print("****************")
        print("Functions info: ")
        print("****************")
        i = 1
        for function in self.get('function'):
            print("")
            print("*** Function: ", i)
            print("")
            function.show()
            i += 1


class BasisSet_C(Structure):
    """
    C interface to the BasisSet class.

    Args:
        basis (BasisSet) : Basis-set of the system.
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
        super(BasisSet_C, self).__init__()

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
                self.l_index[i * 3 + j] = F.get('l')[j]
                self.origin[i * 3 + j] = F.get('origin')[j]

            self.n_prim_cont[i] = F.get('length')
            for j in range(self.n_prim_cont[i]):
                P = F.get('primitive')[j]
                self.exponent[counter] = P.exponent
                self.coefficient[counter] = P.coefficient * P.normalization
                counter += 1
