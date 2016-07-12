# file: timer.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from contextlib import contextmanager
import time


@contextmanager
def timeblock(label):
    start = time.clock()
    try:
        yield
    finally:
        end = time.clock()
        print('\n{0:s} : {1:5.3f} s'.format(label, end - start))
