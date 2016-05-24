# file: test_basis_set.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1

# efposadac@unal.edu.co

import numpy as np
import os
from napmo.system.basis_set import BasisSet


def tests_basis_set_interface():

    test = BasisSet()
    assert test.get('name') == 'user'
    assert test.get('particle') is None
    assert test.get('kind') is None

    basis_file = os.path.join(os.path.dirname(__file__), 'basis.json')
    file = open(basis_file)
    basis_data = file.read().replace('\n', '')
    file.close()

    test.load_gaussian('H', basis_data)
    assert test.get('kind') == 'GTO'
    assert test.get('length') == 1

    value = 0.42377721
    np.testing.assert_allclose(test.compute(), value)

    test = BasisSet('STO-3G')

    assert test.get('name') == 'STO-3G'
    assert test.get('particle') is None
    assert test.get('kind') is None

    test.load_slater('N', basis_data)

    assert test.get('kind') == 'STO'
    assert test.get('length') == 5

    value = np.array([[9.91422306], [-11.93017788], [0.0], [0.0], [0.0]])
    np.testing.assert_allclose(test.compute(), value)

    test2 = BasisSet()
    test2.load_slater('N', basis_data, origin=[1.0, 1.0, 1.0])

    assert test.get('length') == 5
    assert test.get('t_length') == 24

    try:
        print(test2)
    except:
        raise

# tests_basis_set_interface()
