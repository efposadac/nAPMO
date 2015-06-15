# file: atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

import interfaces.elementary_particle
from utilities.databases import AtomicElementsDatabase
from utilities.constants import ANGSTROM_TO_BOHR


class AtomicElement(dict):
    '''
    This class handles all information related to atomic elements Database. (dict)

    Information From: http://physics.nist.gov/constants

    Args:
        symbol (str): Symbol of the atom.
        origin (numpy.ndarray(3)): Origin of the atom (cartesian coordinates)
        BOA (bool, optional): Whether the atom nuclei will be treated in the BOA approach or not. Default is True
        mass_number (int, optional): Mass number of element ``symbol`` :math:`:= A = p^+ + n^o`.
            If 0 the system will choose the most abundant isotope. Default is 0.
        units (str, optional): Units of the origin, valid values are 'Angstroms' or 'Bohr'
    '''
    def __init__(self, symbol, mass_number=0, BOA=True, origin=[0.0, 0.0, 0.0], units='Angstroms'):
        super(AtomicElement, self).__init__()

        assert isinstance(symbol, str)
        assert isinstance(mass_number, int)
        assert isinstance(BOA, bool)
        assert len(origin) == 3
        assert isinstance(units, str)

        try:
            self.update(AtomicElementsDatabase()[symbol])
        except KeyError:
            print('Element: ', symbol, ' not present!')
            raise

        # choose the most abundant mass number (if mass number is not provided)
        aux = AtomicElementsDatabase()['isotopes_' + symbol]['most_abundant']
        self['mass'] = AtomicElementsDatabase()['isotopes_' + symbol][str(aux)]['atomic_weight']

        if mass_number != 0:
            BOA = False

        if not BOA:
            if mass_number != 0:
                try:
                    self.update(AtomicElementsDatabase()['isotopes_' + symbol][str(mass_number)])
                    self['mass'] = self.get('atomic_weight')
                    self['symbol'] = symbol + '_' + str(mass_number)

                except KeyError:
                    print('Mass number: ', str(mass_number), symbol, ' not present!')
                    raise
            else:
                self.update(AtomicElementsDatabase()['isotopes_' + symbol][str(aux)])
                self['symbol'] = symbol + '_' + str(aux)

        # converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'Angstroms':
            origin *= ANGSTROM_TO_BOHR
        self['origin'] = origin

    def isQuantum(self):
        """
        Returns:
            bool: Whether the nuclei of the atomic element is being treated in the BOA framework or not
        """

        output = False

        if 'mass_number' in self:
            output = True

        return output

    def show(self):
        """
        Prints out the information of the object
        """
        print('===================================')
        print('Object: ' + type(self).__name__)
        print('Name: ' + self.get('name'))
        print('Symbol: ' + self.get('symbol'))
        print('Is quantum: ', self.isQuantum())
        print('Z: ', self.get('atomic_number'))
        print('Mass:', self.get('mass'))
        print('origin:', self.get('origin'))
        if 'basis' in self:
            print('Basis set:', self.get('basis').get('name'))
        print('-----------------------------------')
