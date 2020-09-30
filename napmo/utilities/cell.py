# file: cubic_spline.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


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
            self._this = napmo.cext.Cell_new(None, 0)
        else:
            nvec = rvecs.shape[0]
            self._this = napmo.cext.Cell_new(np.PyArray_DATA(rvecs), nvec)

    def __del__(self):
        napmo.cext.Cell_del(self._this)
