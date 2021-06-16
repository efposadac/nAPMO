# file: test_primitive_gaussian.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np

from napmo.gto.primitive_gaussian import PrimitiveGaussian


def test_primitive_gaussian_interface():
    test = PrimitiveGaussian(exponent=0.5)
    assert test.exponent == 0.5
    np.testing.assert_allclose(test.compute(np.zeros([1, 3])), 0.423777208124)
    np.testing.assert_allclose(test.normalization, 0.423777208124)
    np.testing.assert_allclose(test.overlap(test), 1.0)

    test = PrimitiveGaussian(
        exponent=0.5, l=np.array([0, 0, 3], dtype=np.int32))
    np.testing.assert_allclose(test.overlap(test), 1.0)
    np.testing.assert_allclose(test.normalization, 0.309483114995)
    np.testing.assert_allclose(test.compute(np.ones([1, 3])), 0.069055017012)

    try:
        print(test)
    except:
        raise

# test_primitive_gaussian_interface()
