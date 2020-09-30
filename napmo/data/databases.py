# file:  databases.py
# nAPMO package
# Copyright (c) 2014,  Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

import json
import sys
import os


def AtomicElementsDatabase():
    """
    This data base contains all information related to atomic elements.
    Database information taken from: http://physics.nist.gov/constants

    Atomic radii taken from:
    Cordero, et. al. (2008).
    Covalent radii revisited. Dalton Transactions, (21), 2832. http://doi.org/10.1039/b801115j

    .. doctest::

        >>> import napmo
        >>> napmo.AtomicElementsDatabase()['H'] # doctest: +SKIP
        {'atomic_number': 1, 'atomic_radii_2': 0.5858151015155881, \
'electron_affinity': -73, 'atomic_radii': 0.5858151015155881, \
'boiling_point': 20.28, 'charge': -1, 'symbol': 'H', 'density': 0.084, \
'is_quantum': False, 'ionization_energy_1': 1312, 'melting_point': 13.81, \
'electronegativity': 2.1, 'name': 'Hydrogen', 'mass': 1.00794}
    """
    fp = os.path.join(os.path.dirname(__file__), 'atomic_elements.json')

    with open(fp) as f:
        return json.loads(f.read().replace('\n', ''))


def ElementaryParticlesDatabase():
    """
    This database contains all information related to quatum species, ie.
    electrons, muons, etc.

    .. doctest::

        >>> import napmo
        >>> napmo.ElementaryParticlesDatabase()['e-'] # doctest: +SKIP
        {'spin': 0.5, 'name': 'electron', 'category': 'lepton', 'is_quantum': True, \
'mass': 1.0, 'electron': True, 'symbol': 'e-', 'charge': -1.0}
    """
    fp = os.path.join(os.path.dirname(__file__), 'elementary_particles.json')

    with open(fp) as f:
        return json.loads(f.read().replace('\n', ''))


def CouplingConstantsDatabase():
    """
    This database contains all coupling constants for quatum species, ie.
    electrons, muons, etc.

    .. doctest::

        >>> import napmo
        >>> napmo.CouplingConstantsDatabase()['e-'] # doctest: +SKIP
        {'name': 'electron', 'particlesfraction': 0.5, 'kappa': -1.0, \
'eta': 2.0, 'symbol': 'e-', 'lambda': 2.0}
    """
    fp = os.path.join(os.path.dirname(__file__), 'coupling_constants.json')

    with open(fp) as f:
        return json.loads(f.read().replace('\n', ''))
