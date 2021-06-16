# file: test_contracted_gaussian.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np

from napmo.gto.contracted_gaussian import ContractedGaussian


def test_contracted_gaussian_interface():
    test = ContractedGaussian(exponents=[0.5])
    assert test.nprim == 1
    assert test.normalization == 1.0
    print(test)

    np.testing.assert_allclose(np.zeros(3), test.l)
    np.testing.assert_allclose(test.compute(np.zeros([1, 3])), 0.423777208124)
    np.testing.assert_allclose(test.overlap(test), 1.0)

    test = ContractedGaussian(
        exponents=[0.5], l=np.array([0, 0, 3], dtype=np.int32))
    np.testing.assert_allclose(test.overlap(test), 1.0)
    np.testing.assert_allclose(test.compute(np.ones([1, 3])), 0.069055017012)

    try:
        print(test)
    except:
        raise

# test_contracted_gaussian_interface()
