# file: test_elementary_particle.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from napmo.system.elementary_particle import ElementaryParticle
from napmo.gto.basis_set import BasisSet


def test_elementary_particle_interface():
    try:
        a = ElementaryParticle()
        assert False, 'Expecting Failure!'
    except TypeError:
        assert True

    a = ElementaryParticle('e-')
    assert a.get('name') == 'electron'
    assert a.get('symbol') == 'e-'
    assert a.get('charge') == -1.
    assert a.get('mass') == 1.

    try:
        print(a)
    except:
        raise

    try:
        a = ElementaryParticle('MM')
        assert False, 'Expecting Failure!'
    except KeyError:
        assert True

# test_elementary_particle_interface()
