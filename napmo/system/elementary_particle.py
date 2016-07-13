# file: elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


class ElementaryParticle(dict):
    """
    An abstract python interface to create an elementary quantum particle
    i.e Leptons as electron, muon, etc. (dict)

    Database information from:

    From: http://physics.nist.gov/constants

    Args:
        symbol (str): Symbol for the particle.
        origin (ndarray): Origin for the particle (atomic units.)
        units (str, optional): Units of the origin, valid values are
            'ANGSTROMS' or 'BOHR'
    """

    def __init__(self, symbol, origin=np.zeros(3, dtype=np.float64), units='ANGSTROMS'):
        super(ElementaryParticle, self).__init__()

        assert isinstance(symbol, str)

        try:
            self.update(napmo.ElementaryParticlesDatabase()[symbol.lower()])
        except KeyError:
            print('Elementary particle: ', symbol, ' not present!')
            raise

        # converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'ANGSTROMS':
            origin *= napmo.ANGSTROM_TO_BOHR

        self['origin'] = origin

        self.update(napmo.CouplingConstantsDatabase()[symbol.lower()])

    @property
    def is_quantum(self):
        """
        Returns:
            bool: Whether the particle of the atomic element is being treated as punctual particle or not
        """
        return self.get('is_quantum')

    def __repr__(self):

        out = """
==================================================
Object:   {0:9s}
--------------------------------------------------
Name:     {1:9s}
Symbol:   {2:9s}
Size:     {3:<4d}
Category: {4:9s}
Charge:   {5:<5.3f}
Mass:     {6:<5.3f}
Spin:     {7:<5.3f}
Basis:    {8:9s}
""".format(
            type(self).__name__,
            self.get('name'),
            self.get('symbol'),
            self.get('size', 0),
            self.get('category'),
            self.get('charge'),
            self.get('mass'),
            self.get('spin'),
            self.get('basis', {}).get('name', "None")
        )

        if 'origin' in self:
            out += 'Origin:   {0:<5.3f} {1:<5.3f} {2:<5.3f}\n'.format(
                self.get('origin')[0],
                self.get('origin')[1],
                self.get('origin')[2])

        out += '--------------------------------------------------'

        out += ''.join(str(p) for p in self.get('particles', []))

        return out
