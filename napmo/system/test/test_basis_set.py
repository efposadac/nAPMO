# file: test_basis_set.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1

# efposadac@unal.edu.co

import numpy as np
from napmo.system.basis_set import BasisSet, BasisSet_C


def tests_basis_set_interface():

    try:
        test = BasisSet()
        assert False, 'Expecting Failure!'
    except TypeError:
        assert True

    test = BasisSet('STO-3G', 'H')
    assert test.get('name') == 'STO-3G'
    assert test.get('particle') is None
    assert test.get('kind') is None
    assert test.get('length') == 0
    assert test.get('t_length') == 0

    test = BasisSet('STO-3G', 'H', origin=[0.0, 0.0, 0.0])
    assert test.get('length') == 1
    assert test.get('t_length') == 3

    value = np.array([[0.62824688]])
    np.testing.assert_allclose(test.compute(), value)

    test_c = BasisSet_C(test)

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
