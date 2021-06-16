# file: basis_set.py
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
import re
import os


class BasisSet(dict):

    """
    Basis-set interface. (dict)
    """

    def __init__(self, basis_name, particle, origin=None, basis_file=None):
        super(BasisSet, self).__init__()

        self['name'] = basis_name
        self['symbol'] = particle
        self['cont'] = []
        self['nbasis'] = 0

        if origin is not None:

            if not basis_file:
                basis_file = os.path.join(napmo.basis_dir, basis_name)

            basis_data = self.load_file(particle, basis_file)
            self._build(particle, basis_data, origin)

            aux = np.array(
                [cont._this for cont in self.get('cont')], dtype=c_void_p)

            self._this = napmo.cext.BasisSet_new(aux, self.get('nbasis'))

        else:

            self._this = napmo.cext.BasisSet_new_empty()

    def update(self, other):

        self['nbasis'] += other.get('nbasis')
        self['cont'] += other.get('cont')

        napmo.cext.BasisSet_update(self._this, other._this)

    def load_file(self, particle, basis_file):
        """
        Load basis information from file with deMon2k basis-set format. See https://bse.pnl.gov/bse/portal
        for further information.

        Args:
            particle (str): Symbol of the particle.
            basis_file (str): Path to the basis file.

        Return:
            basis_data(list): List of dicts with the basis-set data

            Example::

                [{'angular': 's', 'cont': [0.15432897, 0.15432897, 0.15432897],
                'prim': [3.42525091, 3.42525091, 3.42525091]}]

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

            aux['l'] = int(line[1].strip())
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

    def _build(self, particle, basis_data, origin):
        """
        Load a Gaussian Type Orbital (GTO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            basis_data (list): Information of the basis.
            origin (ndarray) : Center of the functions
        """

        func = [napmo.ContractedGaussian(
            np.array(cont['prim'], dtype=np.float64),
            np.array(cont['cont'], dtype=np.float64),
            np.array(origin, dtype=np.float64),
            np.array([i, k - j, j], dtype=np.int32))
            for cont in basis_data
            for k, i in enumerate(reversed(
                range(cont['l'] + 1)))
            for j in range(k + 1)]

        self['cont'] = func
        self['nbasis'] = len(self.get('cont'))

    def compute(self, coord=np.zeros([1, 3], dtype=np.float64)):
        """
        Compute all the basis-set functions at given ``coord``.

        Args:
            coord (ndarray): coordinates in which the basis set will be
                evaluated. Array shape should be (n, 3)

        """
        n_coord = coord.shape[0]
        output = np.empty([n_coord, self.get('nbasis')])
        napmo.cext.BasisSet_compute(self._this, coord, output, n_coord)

        return output

    def deriv(self, coord=np.zeros([1, 3], dtype=np.float64)):
        """
        Compute first derivative for the basis-set at given ``coord``.

        Args:
            coord (ndarray): coordinates in which the basis set will be
                evaluated. Array shape should be (n, 3)

        """
        n_coord = coord.shape[0]
        output = np.empty([n_coord, self.get('nbasis'), 3])
        napmo.cext.BasisSet_deriv(self._this, coord, output, n_coord)

        return output

    def _lmax(self):
        return np.max(np.array([cont.l.sum() for cont in self.get('cont')]))

    @property
    def lmax(self):
        return self._lmax()

    def __repr__(self):

        out = """
==================================================
Object: {0:9s}       Particle: {1:9s}
--------------------------------------------------
Name: {2:9s} Length: {3:5d}

  {4:<3s} {5:>13s} {6:>13s} {7:>13s}
  ------------------------------------------------
""".format(
            type(self).__name__,
            self.get('symbol', '--'),
            self.get('name', '--'),
            self.get('nbasis', 0),
            "l",
            "zeta",
            "Coeff",
            "Norma"
        )

        out += "\n".join([c._show_compact()
                          for c in self.get('cont', []) if c.l[0] == sum(c.l)])

        out += '--------------------------------------------------'

        return out
