# file: atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np

from napmo.data.databases import AtomicElementsDatabase
from napmo.data.constants import ANGSTROM_TO_BOHR


class AtomicElement(dict):

    '''
    This class handles all information related to atomic elements Database.

    Information From: http://physics.nist.gov/constants

    Args:
        symbol (str): Symbol of the atom.
        origin (ndarray, optional): Origin of the atom (cartesian).
            (Default [0.,0.,0.])
        BOA (bool, optional): Whether the atom nuclei will be treated in the
            BOA framework or not. (Default True)
        mass_number (int, optional): Mass number of element ``symbol``
            :math:`:= A = p^+ + n^o`. If 0 the system will choose the most
            abundant isotope. (Default 0).
        units (str, optional): Units of the origin, valid values are
            'Angstroms' or 'Bohr'
    '''

    def __init__(self, symbol, mass_number=0, BOA=True,
                 origin=np.zeros(3, dtype=np.float64), units='Angstroms'):
        super(AtomicElement, self).__init__()

        assert isinstance(symbol, str)
        assert isinstance(mass_number, int)
        assert isinstance(BOA, bool)
        assert len(origin) == 3
        assert isinstance(units, str)

        data = AtomicElementsDatabase()

        try:
            self.update(data[symbol])
        except KeyError:
            print('Element: ', symbol, ' not present!')
            raise

        # choose the most abundant mass number (if mass number is not provided)
        aux = data['isotopes_' + symbol]['most_abundant']
        self['mass'] = data['isotopes_' + symbol][str(aux)]['atomic_weight']

        if mass_number:
            BOA = False

        if not BOA:
            if mass_number:
                try:
                    self.update(data['isotopes_' + symbol][str(mass_number)])
                    self['mass'] = self.get('atomic_weight')
                    self['symbol'] = symbol + '_' + str(mass_number)

                except KeyError:
                    print('Mass number: ', str(mass_number), symbol,
                          ' not present!')
                    raise
            else:
                self.update(data['isotopes_' + symbol][str(aux)])
                self['symbol'] = symbol + '_' + str(aux)

        # converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'Angstroms':
            origin *= ANGSTROM_TO_BOHR
        self['origin'] = origin

    @property
    def isQuantum(self):
        """
        Returns:
            bool: Whether the nuclei of the atomic element is being treated in
                the BOA framework
        """

        output = False
        if 'mass_number' in self:
            output = True

        return output

    def __repr__(self):
        """
        Prints out the information of the object
        """
        return('===================================\n'
               'Object: ' + type(self).__name__+'\n'
               'Name:   ' + self.get('name')+'\n'
               'Symbol: ' + self.get('symbol')+'\n'
               'Is quantum: ' + str(self.isQuantum)+'\n'
               'Z: ' + str(self.get('atomic_number'))+'\n'
               'Mass:' + str(self.get('mass'))+'\n'
               'origin:' + str(self.get('origin'))+'\n'
               'Basis set:' + self.get('basis', {}).get('name', "None")+'\n'
               '-----------------------------------')
