# file: basis_set.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function
from ctypes import *

import numpy as np
import napmo as nap
import re
import os


class BasisSet(dict):

    """
    Basis-set interface. (dict)

    Basis-set manager
    """

    def __init__(self, basis_name, particle, origin=None, basis_file=None):
        super(BasisSet, self).__init__()
        self['name'] = basis_name
        self['symbol'] = particle
        self['function'] = []
        self['length'] = 0
        self['t_length'] = 0

        if origin is not None:
            if not basis_file:
                basis_file = os.path.join(nap.basis_dir, basis_name)

            basis_data = self.load_file(particle, basis_file)
            self.load_gaussian(particle, basis_data, origin)

    def update(self, other):
        self['function'] += other.get('function', [])
        self['length'] += other.get('length', 0)
        self['t_length'] += other.get('t_length', 0)

    def load_file(self, particle, basis_file):
        """
        Load basis information from file with deMon2k basis-set format. See https://bse.pnl.gov/bse/portal
        for further information.

        Args:
            particle (str): Symbol of the particle.
            basis_file (str): Path to the basis file.

        Return:
            basis_data(list): dict list with the basis-set data

            Example:

            [{'angular': 's', 'cont': [0.15432897, 0.15432897, 0.15432897], 'prim': [3.42525091, 3.42525091, 3.42525091]}]

        """

        basis_data = []

        with open(basis_file, 'r') as f:
            data = f.read()

        # Clean commentaries and blank lines
        data = re.sub('#.*\n?', '', data)
        data = re.sub('(?imu)^\s*\n', '', data)

        # Find atom
        particle = particle.replace('+', '\+')
        particle = particle.replace('-', '\-')
        pos = re.search('O-.*\s' + particle + '\s', data)

        if pos is None:
            try:
                raise ValueError('Basis for particle ' +
                                 particle + ' not found!')
            except:
                print('Check you basis file: ', basis_file)
                raise

        data = data[pos.start():].splitlines()

        line_index = 2
        for cont in range(int(data[1].strip())):
            aux = {}
            line = data[line_index].strip().split()

            aux['angular'] = int(line[1].strip())
            aux['cont'] = []
            aux['prim'] = []

            line_index += 1

            for prim in range(int(line[2].strip())):
                line = data[line_index].strip().split()

                aux['prim'].append(float(line[0].strip()))
                aux['cont'].append(float(line[1].strip()))

                line_index += 1

            basis_data.append(aux)

        return basis_data

    def load_gaussian(self, particle, basis_data, origin):
        """
        Load a Gaussian Type Orbital (GTO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            basis_data (list): Information of the basis.
        """

        func = [nap.ContractedGaussian(
            np.array(cont['prim'], dtype=np.float64),
            np.array(cont['cont'], dtype=np.float64),
            np.array(origin, dtype=np.float64),
            np.array([i, k - j, j], dtype=np.int32))
            for cont in basis_data
            for k, i in enumerate(reversed(
                range(cont['angular'] + 1)))
            for j in range(k + 1)]

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
        output = (func.compute(coord) for func in self.get('function'))
        return np.column_stack(output)

    def __repr__(self):

        out = """
==================================================
Object: {0:9s}       Particle: {1:9s}
--------------------------------------------------

Name: {2:9s} Length: {3:5d}

  {4:<3s} {5:>10s} {6:>10s} {7:>10s}
  ------------------------------------
""".format(
            type(self).__name__,
            self.get('symbol', '--'),
            self.get('name', '--'),
            self.get('length', 0),
            "l",
            "zeta",
            "Coeff",
            "Norma"
        )

        out += "\n".join([p._show_compact()
                          for p in self.get('function', []) if p.l[0] == sum(p.l)])

        out += '--------------------------------------------------'

        return out


class BasisSet_C(Structure):

    """
    C interface to the BasisSet class.

    Args:
        basis(BasisSet): Basis set of the system.
    """
    _fields_ = [
        ("n_cont", c_int),
        ("n_prim_cont", POINTER(c_int)),
        ("prim_index", POINTER(c_int)),
        ("l_index", POINTER(c_int)),
        ("origin", POINTER(c_double)),
        ("normalization", POINTER(c_double)),
        ("exponent", POINTER(c_double)),
        ("coefficient", POINTER(c_double)),
        ("p_normalization", POINTER(c_double)),
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

        size = (c_double * basis.get('t_length'))()
        self.p_normalization = cast(size, POINTER(c_double))

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
                self.coefficient[counter] = P.coefficient
                self.p_normalization[counter] = P.normalization
                counter += 1
