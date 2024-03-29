# file: auxiliary_basis.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np
from numpy import linalg as LA

import napmo


class AuxiliaryBasis(object):
    """
    This class handles the auxiliary basis for grid-based HF

    Args:
        nspecies (int): Total number of unique species in the system

    """

    def __init__(self, grid):
        super(AuxiliaryBasis, self).__init__()

        assert isinstance(grid, napmo.BeckeGrid)

        self._grid = grid
        self._ncenter = grid.ncenter
        self._this = napmo.cext.AuxiliaryBasis_new(self.grid._this)
        self._nao = napmo.cext.AuxiliaryBasis_get_nao(self._this)

    @property
    def basis(self):
        """
        Returns the auxiliary basis
        """
        ptr = napmo.cext.AuxiliaryBasis_get_aobasis(self._this)
        return np.ctypeslib.as_array(ptr, shape=(self.grid.size, self.nao)).T

    @property
    def ncenter(self):
        return self._ncenter

    @property
    def nao(self):
        return self._nao

    @property
    def grid(self):
        return self._grid
