# file: test_basis_set.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0

# efposadac@sissa.it

import os
import sys
import numpy as np

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.basis_set import BasisSet


def tests_basis_set_interface():

    test = BasisSet()
    assert test.get('name') == 'user'
    assert test.get('particle') == None
    assert test.get('kind') == None

    file = open('basis.json')
    basis_data = file.read().replace('\n', '')
    file.close()

    test.load_gaussian('H', basis_data)
    assert test.get('kind') == 'GTO'
    assert test.get('length') == 1

    try:
        test.show()
    except:
        raise

    value = 0.42377721
    np.testing.assert_allclose(test.compute(), value)

    test = BasisSet('STO-3G')

    assert test.get('name') == 'STO-3G'
    assert test.get('particle') == None
    assert test.get('kind') == None

    test.load_slater('N', basis_data)

    assert test.get('kind') == 'STO'
    assert test.get('length') == 5

    value = np.array([9.91422306, -11.93017788, 0.0, 0.0, 0.0])
    np.testing.assert_allclose(test.compute(), value)

    test2 = BasisSet()
    test2.load_slater('N', basis_data, origin=[1.0, 1.0, 1.0])

    test += test2
    assert test.get('length') == 10
    assert test.get('t_length') == 48

    try:
        test2.show()
    except:
        raise

    try:
        test2.show_json()
    except:
        raise

# tests_basis_set_interface()
