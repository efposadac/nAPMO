# file: test_timer.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

import napmo
import numpy as np


def test_timer():
    with napmo.runtime.timeblock('Test'):
        for i in range(100000):
            x = np.sqrt(i**1.5)

    napmo.runtime.show_block('Test')

    with napmo.runtime.timeblock('Test'):
        for i in range(100000):
            x = np.sqrt(i**1.5)

    napmo.runtime.show_block('Test')

    with napmo.runtime.timeblock('Test2'):
        for i in range(100000):
            x = np.sqrt(i**1.5)

    napmo.runtime.show_block('Test2')

    napmo.runtime.show_summary()

# test_timer()
