# file: test_atomic_element.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import numpy as np

from napmo.system.atomic_element import AtomicElement
from napmo.gto.basis_set import BasisSet


def test_atomic_element_interface():

    try:
        a = AtomicElement()
        assert False, 'Expecting Failure!'
    except TypeError:
        assert True

    try:
        a = AtomicElement('UNKNOWN')
        assert False, 'Expecting Failure!'
    except KeyError:
        assert True

    a = AtomicElement('C')

    assert a.get('name') == 'Carbon'
    assert a.is_quantum is False, 'Expecting Failure!'

    for i in range(3):
        assert a.get('origin')[i] == 0.

    a = AtomicElement('H_1')
    assert a.get('mass_number') == 1
    assert a.is_quantum is True

    try:
        a = AtomicElement('H_100')
        assert False, 'Expecting Failure!'
    except KeyError:
        assert True

    try:
        print(a)
    except:
        raise

# test_atomic_element_interface()
