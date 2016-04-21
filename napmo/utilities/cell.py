# file: cubic_spline.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import numpy.ctypeslib as npct

from napmo.system.cext import napmo_library as nl


class Cell(object):

    """
    Representation of periodic boundary conditions.

       0, 1, 2 and 3 dimensional systems are supported. The cell vectors don't
       need to be orthogonal.
    """

    def __init__(self, rvecs=None, initvoid=False):
        super(Cell, self).__init__()
        if initvoid:
            self._this = None
        elif rvecs is None:
            self._this = nl.Cell_new(None, 0)
        else:
            nvec = rvecs.shape[0]
            self._this = nl.Cell_new(np.PyArray_DATA(rvecs), nvec)

    def __del__(self):
        nl.Cell_del(self._this)

nl.Cell_new.restype = c_void_p
nl.Cell_new.argtypes = [c_void_p, c_int]

nl.Cell_del.restype = None
nl.Cell_del.argtypes = [c_void_p]
