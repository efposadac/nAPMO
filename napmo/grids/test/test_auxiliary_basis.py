# file: test_auxiliary_basis.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA
from scipy.linalg import schur, eigvals
import re


def test_auxiliary_basis():
    molecule = napmo.MolecularSystem()
    molecule.add_atom(
        "He", [0.0, 0.0, 0.0], basis_name='STO-3G'
    )

    angularPoints = 110
    radialPoints = 10
    lmax = int(napmo.lebedev_get_order(angularPoints) / 2.0)
    grid = napmo.BeckeGrid(molecule.get('e-'), radialPoints, angularPoints)

    with napmo.runtime.timeblock('C'):
        aobasis = napmo.AuxiliaryBasis(grid)

    assert(aobasis.ncenter == 1)
    assert(aobasis.nao == 9)
    assert(aobasis.basis.shape == (9, grid.size))


if __name__ == '__main__':
    test_auxiliary_basis()
