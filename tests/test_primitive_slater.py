# file: test_primitive_slater.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np

from napmo.interfaces.primitive_slater import PrimitiveSlater


def test_primitive_slater_interface():
    test = PrimitiveSlater(exponent=0.5)
    assert test.get('exponent') == 0.5
    np.testing.assert_allclose(test.compute(np.zeros(3)), 0.19947114020071635)
    np.testing.assert_allclose(test.get('normalization'), 0.7071067811865476)
    np.testing.assert_allclose(test.overlap(test), 1.0)

    test = PrimitiveSlater(exponent=0.5, l=3)
    np.testing.assert_allclose(test.overlap(test), 1.0)
    np.testing.assert_allclose(test.get('normalization'), 0.7071067811865476)
    np.testing.assert_allclose(test.compute(np.ones(3)), -0.03924135968420833)

test_primitive_slater_interface()
