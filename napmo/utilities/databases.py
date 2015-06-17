# file:  databases.py
# nAPMO package
# Copyright (c) 2014,  Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import json
import sys
import os


def AtomicElementsDatabase():
    pass
    """
    AtomicElementsDatabase
    ======================

    This data base contains all information related to atomic elements.
    Database information taken from: http://physics.nist.gov/constants

    Atomic radii from J. C. Slater, J. Chem. Phys. 41, 3199 (1964).

    Example:

        >>> import utilities.databases
        >>> utilities.databases.AtomicElementsDatabase()['H'] # doctest: +SKIP
        {'vanderwaals_radius': 1.2, 'atomic_number': 1, 'symbol': 'H', 'covalent_radius': 0.3, \
    'atomic_radii': 0.35, 'ionization_energy_1': 1312, 'melting_point': 13.81, 'electron_affinity': -73, \
    'electronegativity': 2.1, 'name': 'Hydrogen', 'density': 0.084, 'boiling_point': 20.28, 'mass': 1.00794}

    """
    fp = os.path.join(os.path.dirname(__file__), 'atomic_elements.json')
    file = open(fp)
    data = file.read().replace('\n', '')
    file.close()

    return json.loads(data)


def ElementaryParticlesDatabase():
    """
        ElementaryParticlesDatabase
    ===========================

    This database contains all information related to quatum species, ie. electrons, muons, etc.

    Example:

        >>> import utilities.databases
        >>> utilities.databases.ElementaryParticlesDatabase()['e-'] # doctest: +SKIP
        {'category': 'lepton', 'name': 'electron', 'symbol': 'e-', 'charge': -1.0, 'mass': 1.0, 'spin': 0.5}

    UnitsDatabase
    =============

    This database contains factors to covert between Bohr and Angstroms.

    Example:

        >>> import utilities.databases
        >>> utilities.databases.UnitsDatabase['Angstroms']
        0.52917721092
    """
    fp = os.path.join(os.path.dirname(__file__), 'elementary_particles.json')
    file = open(fp)
    data = file.read().replace('\n', '')
    file.close()

    return json.loads(data)
