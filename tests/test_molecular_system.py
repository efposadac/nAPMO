# file: test_molecular_system.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

import os
import sys

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.molecular_system import *


def test_molecular_system_interface():

    a = MolecularSystem()
    a.add_atom('C', [0., 0., 200.])

    assert a.n_elementary_particles() == 1  # e-
    assert a.n_atoms() == 1  # C

    a.add_elementary_particle('e-', [0, 0, 100])
    assert a.n_elementary_particles() == 1
    assert a.n_particles('e-') == 7

    a.add_elementary_particle('u-', [0, 0, 0])

    assert a.n_elementary_particles() == 2

    try:
        a.n_particles('x-')
        assert False, 'Failure expected!'
    except:
        pass

    try:
        a.show()
        assert True
    except:
        pass

    try:
        test = a['e+']
        assert False, 'Failure expected!'
    except:
        pass

    assert a['u-'].peek().get('name') == 'muon'

    assert a['atoms'][0].isQuantum() is False

    a.add_atom('H', [0., 0., 0.], BOA=False)

    assert a['atoms'][1].isQuantum() is True

test_molecular_system_interface()
