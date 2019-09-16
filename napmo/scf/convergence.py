# file: convergence.py
# nAPMO package
# Copyright (c) 2017, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo

from scipy import linalg as SLA


class Convergence(object):
    """
    Implements several converge algorithms
    """
    def __init__(self, F, D):
        super(Convergence, self).__init__()

        self._ndim = F.shape[0]

        self.IF = np.zeros([self._ndim, self._ndim])
        self.ID = np.zeros([self._ndim, self._ndim])
        self.X = np.zeros([self._ndim, self._ndim])

        self.IF[:] = F
        self.ID[:] = D

    def damping(self, NF, ND):
        """
        Damping implementation
        """

        # Computes the effect of the changes in the Density matrix
        self.X[:] = self.IF
        self.X[:] = self.X.dot(ND - self.ID)

        dens_effect = -0.5 * self.X.trace()

        # Computes the effect of the changes in the Fock and Density matrices
        self.X[:] = (NF - self.IF).dot(ND - self.ID)

        dens_fock_effect = self.X.trace()

        # Check whether to change the Fock matrix or not
        output = NF
        damping_factor = 0.0

        if dens_fock_effect <= dens_effect:
            self.IF[:] = NF
            self.ID[:] = ND
        else:
            damping_factor = dens_effect / dens_fock_effect
            self.IF[:] += (damping_factor * (NF - self.IF))
            output = self.IF

            self.ID[:] += damping_factor * (ND - self.ID)

        return output

    @property
    def ndim(self):
        return self._ndim
