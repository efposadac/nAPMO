# file: c_binding.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

from ctypes import *
import os

try:
    lp = os.path.join(os.path.dirname(__file__), 'libnapmo.so')
    napmo_library = CDLL(lp)
except OSError:
    lp = os.path.join(os.path.dirname(__file__), '../utilities/libnapmo.so')
    napmo_library = CDLL(lp)
