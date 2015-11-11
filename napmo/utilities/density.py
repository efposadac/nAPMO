# file: density.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy.ctypeslib as npct
import numpy as np
from ctypes import *

from napmo.system.cext import napmo_library
from napmo.system.basis_set import BasisSet_C

array_1d_double = npct.ndpointer(
    dtype=np.double, ndim=1, flags='CONTIGUOUS')

array_2d_double = npct.ndpointer(
    dtype=np.double, ndim=2, flags='CONTIGUOUS')

napmo_library.density_gto.restype = None
napmo_library.density_gto.argtypes = [
    POINTER(BasisSet_C),
    array_2d_double,
    array_2d_double,
    array_1d_double,
    c_int
]


def density_full_from_matrix_gto(density_file, basis, coords):
    """
    Calculates the value of density :math:`\\rho` at point ``coord``


    Args:
        density_file (str): name of the file which contains the density matrix.
        basis (BasisSet_C): Basis set allocated in C.
        coords (ndarray): array with the points to calculate the density. ``ndim=2, shape=(n,3)``

    Returns:
        rho (ndarray): Array with the value of density :math:`\\rho` at ``coords`` points.
    """
    size = coords.shape[0]
    P = np.array(np.loadtxt(density_file), order='C', dtype=np.float64)
    rho = np.empty(size, dtype=np.float64)
    napmo_library.density_gto(byref(basis), coords, P, rho, size)
    return rho
