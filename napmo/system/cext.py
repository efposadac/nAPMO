# file: cext.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from ctypes import *
import os

try:
    lp = os.path.join(os.path.dirname(__file__), 'libnapmo.so')
    napmo_library = CDLL(lp)
except OSError:
    lp = os.path.join(os.path.dirname(__file__), '../utilities/libnapmo.so')
    napmo_library = CDLL(lp)
