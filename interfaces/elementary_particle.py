# file: elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from utilities.databases import ElementaryParticlesDatabase


class ElementaryParticle(dict):
    """
    An abstract python interface to create an elementary quantum particle
    i.e Leptons as electron, muon, etc. (dict)

    Database information from:

    From: http://physics.nist.gov/constants

    Args:
            symbol (str): Symbol for the particle.
            origin (numpy.ndarray(3)): Origin for the particle (atomic units.)
    """

    def __init__(self, symbol='null', origin=[0.0, 0.0, 0.0]):
        super(ElementaryParticle, self).__init__()

        assert isinstance(symbol, str)
        assert len(origin) == 3

        try:
            self.update(ElementaryParticlesDatabase()[symbol])
        except KeyError:
            print('Elementary particle: ', symbol, ' not present!, creating one.')
            self.update(ElementaryParticlesDatabase()['user'])
            self['symbol'] = symbol

        self['origin'] = np.array(origin, dtype=np.float64)

    def show(self):
        """
        Shows the information of the objects
        """
        print('===================================')
        print('Object: ' + type(self).__name__)
        print('Name: ' + self.get('name'))
        print('Symbol: ' + self.get('symbol'))
        print('Category: ' + self.get('category'))
        print('Charge:', self.get('charge'))
        print('Mass:', self.get('mass'))
        print('Spin:', self.get('spin'))
        print('origin:', self.get('origin'))
        if 'basis' in self:
            print('Basis set:', self.get('basis')['name'])
        print('-----------------------------------')
