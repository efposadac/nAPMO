# file: test_primitive_gaussian.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np

from napmo.interfaces.primitive_gaussian import PrimitiveGaussian


def test_primitive_gaussian_interface():
    test = PrimitiveGaussian(exponent=0.5)
    assert test.get('exponent') == 0.5
    np.testing.assert_allclose(test.compute(np.zeros(3)), 0.423777208124)
    np.testing.assert_allclose(test.get('normalization'), 0.423777208124)
    np.testing.assert_allclose(test.overlap(test), 1.0)

    test = PrimitiveGaussian(exponent=0.5, l=[0, 0, 3])
    np.testing.assert_allclose(test.overlap(test), 1.0)
    np.testing.assert_allclose(test.get('normalization'), 0.309483114995)
    np.testing.assert_allclose(test.compute(np.ones(3)), 0.069055017012)

test_primitive_gaussian_interface()
