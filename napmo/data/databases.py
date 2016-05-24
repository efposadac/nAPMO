# file:  databases.py
# nAPMO package
# Copyright (c) 2014,  Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import json
import sys
import os


def AtomicElementsDatabase():
    """
    This data base contains all information related to atomic elements.
    Database information taken from: http://physics.nist.gov/constants

    Atomic radii from J. C. Slater, J. Chem. Phys. 41, 3199 (1964).

    """
    fp = os.path.join(os.path.dirname(__file__), 'atomic_elements.json')

    with open(fp) as f:
        return json.loads(f.read().replace('\n', ''))


def ElementaryParticlesDatabase():
    """
    This database contains all information related to quatum species, ie.
    electrons, muons, etc.
    """
    fp = os.path.join(os.path.dirname(__file__), 'elementary_particles.json')

    with open(fp) as f:
        return json.loads(f.read().replace('\n', ''))
