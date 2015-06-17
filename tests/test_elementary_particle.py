# file: test_elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from napmo.interfaces.elementary_particle import ElementaryParticle
from napmo.interfaces.basis_set import BasisSet


def test_elementary_particle_interface():
    a = ElementaryParticle()
    assert a.get('name') == 'user'

    a['name'] = 'electron'
    assert a.get('name') == 'electron'

    a = ElementaryParticle('e-', origin=[0, 1, 2])
    assert a.get('name') == 'electron'
    assert a.get('symbol') == 'e-'
    assert a.get('charge') == -1.
    assert a.get('mass') == 1.

    for i in range(3):
        assert a.get('origin')[i] == float(i)

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
