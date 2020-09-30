# file: test_lebedev.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

import napmo.grids.lebedev as nint
import numpy as np


def test_lebedev_get_order():
    assert nint.lebedev_get_order(110) == 17

    try:
        a = nint.lebedev_get_order(2)
        assert False, 'Expecting Failure!'
    except:
        assert True

# test_lebedev_get_order()
