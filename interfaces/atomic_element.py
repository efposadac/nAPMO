# file: atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from utilities.databases import AtomicElementsDatabase
from utilities.databases import UnitsDatabase
from copy import deepcopy


class AtomicElement(object):
    '''This class handles all information related to atomic elements
    Database information From: http://physics.nist.gov/constants
    symbol : Atomic element symbol
    mass_number : Mass number of element 'symbol' := A = p+ + n
    BOA : treat nuclei in the Born-Oppenheimer approximation framework.
    '''
    def __init__(self, symbol, mass_number=0, BOA=True, position=[0.0, 0.0, 0.0], units='Angstroms'):
        '''Returns an object with all information of the atomic element symbol.
        data: python dictionary with information related to atom 'symbol'
        '''
        super(AtomicElement, self).__init__()

        assert isinstance(symbol, str)
        assert isinstance(mass_number, int)
        assert isinstance(BOA, bool)
        assert len(position) == 3
        assert isinstance(units, str)

        try:
            self.data = deepcopy(AtomicElementsDatabase[symbol])
        except KeyError:
            print 'Element: ', symbol, ' not present!'
            raise

        # most abundant mass number (if mass number is not provided)
        aux = AtomicElementsDatabase['isotopes_'+symbol]['most_abundant']
        self.set('mass', AtomicElementsDatabase['isotopes_'+symbol][aux]['atomicWeight'])

        if mass_number != 0:
            BOA = False

        if not BOA:
            if mass_number != 0:
                try:
                    self.data.update(AtomicElementsDatabase['isotopes_'+symbol][mass_number])
                    self.set('mass', self.data['atomicWeight'])
                    self.set('symbol', symbol+'_'+str(mass_number))

                except KeyError:
                    print 'Mass number: ', str(mass_number), symbol, ' not present!'
                    raise
            else:
                self.data.update(AtomicElementsDatabase['isotopes_'+symbol][aux])
                self.set('symbol', symbol+'_'+str(aux))

        # converting to Bohr
        position = np.array(position, dtype=np.float64)
        if units == 'Angstroms':
            position *= UnitsDatabase['Bohr']
        self.set('position', position)

    def isQuantum(self):
        """ Returns whether the nuclei of the atomic element is being
        treated in the BOA framework or not"""
        try:
            self.data['mass_number']
            return True
        except KeyError:
            return False

    def get(self, key):
        """Returns the value stored in key
        """
        assert isinstance(key, str)

        try:
            return self.data[key]
        except KeyError:
            raise

    def set(self, key, value):
        """set a new value.
        """
        # TODO: Improve this function!
        self.data[key] = value

    def show(self):
        """Shows the information of the object
        """
        print '==================================='
        print 'Object: '+type(self).__name__
        print 'Name: '+self.get('name')
        print 'Symbol: '+self.get('symbol')
        print 'Is quantum: ', self.isQuantum()
        print 'Z: ', self.get('atomicNumber')
        print 'Mass:', self.get('mass')
        print 'Position:', self.get('position')
        print '-----------------------------------'
