# file: test_atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

import numpy as np

import os
import sys

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.atomic_element import AtomicElement


def test_atomic_element_interface():

    try:
        a = AtomicElement()
        assert False, 'Failure expected!'
    except:
        pass

    try:
        a = AtomicElement('UNKNOWN')
        assert False, 'Failure expected'
    except:
        pass

    a = AtomicElement('C')

    assert a.get('name') == 'Carbon'
    assert a.isQuantum() is False

    for i in range(3):
        assert a.get('origin')[i] == 0.

    a = AtomicElement('H', BOA=False)

    assert a.isQuantum() is True
    assert a.get('mass_number') == 1

    a['mass_number'] = 2
    assert a.get('mass_number') == 2

    a['symbol'] = 'none'
    assert a.get('symbol') == 'none'

    a = AtomicElement('H', mass_number=3)
    assert a.isQuantum() is True

    try:
        a = AtomicElement('H', mass_number=100)
        assert False, 'Failure expected'
    except:
        pass

    try:
        a.get('address')
        assert False, 'Failure expected'
    except:
        pass

    a['address'] = 'sky'
    assert a.get('address') == 'sky'

    try:
        a.show()
        assert True
    except:
        pass

# test_atomic_element_interface()
