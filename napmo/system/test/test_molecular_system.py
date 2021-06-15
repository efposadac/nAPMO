# file: test_molecular_system.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from napmo.system.molecular_system import MolecularSystem


def test_molecular_system_interface():

    a = MolecularSystem()
    a.add_atom('C', [0., 0., 200.])
    assert a.get('e-').get('particles')[0].is_quantum is False

    assert a.size_species == 1  # e-

    a.add_elementary_particle('e-', [0, 0, 100])
    assert a.get('e-').get('particles')[1].is_quantum is False
    assert a.size_species == 1
    assert a.size_particles('e-') == 7

    a.add_elementary_particle('u-', [0, 0, 0])

    assert a.size_species == 2
    assert a.size_particles('x-') == 0

    a.add_nucleus('H_2', [0, 0, 0])
    assert a.size_species == 3

    a.add_nucleus('H', [0, 0, 0])
    assert a.get('e-').get('particles')[2].is_quantum is False
    assert a.size_species == 3

    try:
        test = a['e+']
        assert False, 'Expecting Failure!'
    except KeyError:
        assert True

    assert a.get('u-').get('name') == 'muon'

    a.add_atom('He', [0., 0., 0.], quantum=True)

    assert a.get('e-').get('particles')[3].is_quantum is True

    print(a.keys())

    try:
        print(a)
    except:
        raise

    # More real test
    particle = "N"
    basis_name = "STO-3G"

    molecule = MolecularSystem()

    molecule.add_atom(
        particle, [0.000000, 0.000000, 0.70997005], basis_name=basis_name)

    molecule.add_atom(
        particle, [0.000000, 0.000000, -0.70997005], basis_name=basis_name)

    assert molecule.size_particles('e-') == 14
    assert molecule.size_species == 1

    charges = {'e-': {'charge': 1, 'multiplicity': 1}}
    try:
        molecule.set_charges(charges, open_shell=True)
        assert False, 'Expecting Failure'
    except ValueError:
        assert True

    print(molecule)

    particle = "N"
    basis_name = "STO-3G"

    molecule = MolecularSystem()

    molecule.add_atom(
        particle, [0.000000, 0.000000, 0.70997005], basis_name=basis_name, quantum=True)

    molecule.add_atom(
        particle, [0.000000, 0.000000, -0.70997005], basis_name=basis_name, quantum=True)

    molecule.add_elementary_particle(
        'U+', [0.000000, 0.000000, -0.50997005], basis_name='NAKAI-5-SP')

    molecule.add_nucleus(
        'H_1', [0.000000, 0.000000, -0.70997005], basis_name='DZSPDN')

    molecule.add_nucleus(
        'H_1', [0.000000, 0.000000, 0.70997005], basis_name='DZSPDN')

    assert molecule.size_particles('e-') == 14
    assert molecule.size_species == 3

    basis = molecule.get_basis('H_1')
    basis = molecule.get_basis('N')

    charges = {'e-': {'charge': 1, 'multiplicity': 2},
               'U+': {'charge': 1}}

    molecule.set_charges(charges, open_shell=True)

    charges = {'e-': {'charge': 1, 'multiplicity': 2}}
    molecule.set_charges(charges, open_shell=True)

    print(molecule)

    particle = "N"
    basis_name = "STO-3G"

    molecule = MolecularSystem()

    molecule.add_atom(
        particle, [0.000000, 0.000000, 0.70997005], basis_name=basis_name, quantum=True)

    charges = {'e-': {'charge': 1, 'multiplicity': 1}}
    molecule.set_charges(charges, open_shell=False)


# test_molecular_system_interface()
