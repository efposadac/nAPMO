# file: nwavefunction.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


def print_matrix(m):
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            print('{0:>7.4f} '.format(m[i, j]), end=' ')
        print('')


class NWaveFunction(object):
    """
    Numerical WaveFunction class

    This class is initiated with some matrices obtained from a minimal basis analytical
    calculation.

    Args:
        psi (WaveFunction) : Optimized wavefunction object
    """

    def __init__(self, psi, nrad, nang):
        super(NWaveFunction, self).__init__()
        self._apsi = psi
        self._mgrid = napmo.BeckeGrid(psi.species, nrad, nang)

        self._mgrid.show()

        # Expand density
        print_matrix(psi.D)

        print(psi.species.get('basis').compute(self._mgrid.points))

        print('here we go')
