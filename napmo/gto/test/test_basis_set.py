# file: test_basis_set.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0

# fernando.posada@temple.edu

import numpy as np
from napmo.gto.basis_set import BasisSet


def tests_basis_set_interface():

    try:
        test = BasisSet()
        assert False, 'Expecting Failure!'
    except TypeError:
        assert True

    test = BasisSet('STO-3G', 'H')
    assert test.get('name') == 'STO-3G'
    assert test.get('particle') is None
    assert test.get('nbasis') == 0

    test = BasisSet('STO-3G', 'H', origin=[0.0, 0.0, 0.0])
    assert test.get('nbasis') == 1

    value = np.array([[0.62824688]])
    np.testing.assert_allclose(test.compute(), value)

    try:
        print(test)
    except:
        raise

    try:
        test = BasisSet('STO-3G', 'HH', origin=[0.0, 0.0, 0.0])
        assert False, 'Expecting Failure!'
    except ValueError:
        assert True

# tests_basis_set_interface()
