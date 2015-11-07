# file: test_angular_quadratures.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import napmo.utilities.angular_quadratures as nint
import numpy as np


def test_angular_quadratures_lebedev_get_order():
    assert nint.lebedev_get_order(110) == 17

    try:
        a = nint.lebedev_get_order(2)
        assert False, 'Expecting Failure!'
    except:
        assert True

# test_angular_quadratures_lebedev_get_order()
