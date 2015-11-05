# file: test_molecular_system.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
import os

from napmo.system.molecular_system import MolecularSystem


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
        assert False, 'Expecting Failure!'
    except:
        assert True

    try:
        a.show()
    except:
        raise

    try:
        test = a['e+']
        assert False, 'Expecting Failure!'
    except:
        assert True

    assert a.get('u-')[-1].get('name') == 'muon'

    assert a.get('atoms')[0].isQuantum() is False

    a.add_atom('H', [0., 0., 0.], BOA=False)

    assert a.get('atoms')[1].isQuantum() is True

    # More real test
    particle = "N"
    basis_name = "STO-3G"
    basis_file = os.path.join(
        os.path.dirname(__file__), "basis.json")
    basis_kind = "GTO"

    molecule = MolecularSystem()

    molecule.add_atom(
        particle, [0.000000, 0.000000, 0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    molecule.add_atom(
        particle, [0.000000, 0.000000, -0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    molecule.add_elementary_particle(
        'u+', [0.000000, 0.000000, -0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    assert molecule.n_particles('e-') == 14
    assert molecule.n_elementary_particles() == 2
    assert molecule.n_atoms() == 2
    molecule.show()

    particle = "N"
    basis_name = "STO-3G"
    basis_file = basis_file = os.path.join(
        os.path.dirname(__file__), "basis.json")
    basis_kind = "STO"

    molecule = MolecularSystem()

    molecule.add_atom(
        particle, [0.000000, 0.000000, 0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    molecule.add_atom(
        particle, [0.000000, 0.000000, -0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    molecule.add_elementary_particle(
        'u-', [0.000000, 0.000000, -0.70997005],
        basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
    )

    assert molecule.n_particles('e-') == 14
    assert molecule.n_particles('u-') == 1
    assert molecule.n_elementary_particles() == 2
    assert molecule.n_atoms() == 2

    try:
        molecule.show()
    except:
        raise

# test_molecular_system_interface()
