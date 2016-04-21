# file: test_elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from napmo.system.elementary_particle import ElementaryParticle
from napmo.system.basis_set import BasisSet


def test_elementary_particle_interface():
    a = ElementaryParticle()
    assert a.get('name') == 'user'

    a['name'] = 'electron'
    assert a.get('name') == 'electron'

    a = ElementaryParticle('e-')
    assert a.get('name') == 'electron'
    assert a.get('symbol') == 'e-'
    assert a.get('charge') == -1.
    assert a.get('mass') == 1.

    try:
        a.get('address')
        assert False, 'Expecting Failure!'
    except:
        assert True

    basis = BasisSet()
    a['basis'] = basis

    try:
        a.show()
    except:
        raise

# test_elementary_particle_interface()
