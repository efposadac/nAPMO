# file: test_atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import numpy as np

from napmo.system.atomic_element import AtomicElement
from napmo.system.basis_set import BasisSet


def test_atomic_element_interface():

    try:
        a = AtomicElement()
        assert False, 'Expecting Failure!'
    except:
        assert True

    try:
        a = AtomicElement('UNKNOWN')
        assert False, 'Expecting Failure!'
    except:
        assert True

    a = AtomicElement('C')

    assert a.get('name') == 'Carbon'
    assert a.isQuantum is False, 'Expecting Failure!'

    for i in range(3):
        assert a.get('origin')[i] == 0.

    a = AtomicElement('H', BOA=False)

    assert a.isQuantum is True
    assert a.get('mass_number') == 1

    a['mass_number'] = 2
    assert a.get('mass_number') == 2

    a['symbol'] = 'none'
    assert a.get('symbol') == 'none'

    a = AtomicElement('H', mass_number=3)
    assert a.isQuantum is True

    try:
        a = AtomicElement('H', mass_number=100)
        assert False, 'Expecting Failure!'
    except:
        assert True

    try:
        a.get('address')
        assert False, 'Expecting Failure!'
    except:
        assert True

    a['address'] = 'sky'
    assert a.get('address') == 'sky'

    basis = BasisSet()
    a['basis'] = basis

# test_atomic_element_interface()
