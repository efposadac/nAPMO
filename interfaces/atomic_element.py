# file: atomic_element.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from databases import AtomicElementsDatabase
from databases import UnitsDatabase


class AtomicElement(object):
    '''This class handles all information related to atomic elements
    
    Database information From: http://physics.nist.gov/constants

    symbol : Atomic element symbol
    mass_number : Mass number of element 'symbol' := A = p+ + n
    BOA : treat nuclei in the Born-Oppenheimer approximation framework.

    '''
    def __init__(self, symbol, mass_number=0, BOA=True, position=[0.0,0.0,0.0], units='Angstroms'):
        '''Returns an object with all information of the atomic element symbol.
        data: python dictionary with information related to atom 'symbol'
        '''
        super(AtomicElement, self).__init__()

        try:
            self.data = AtomicElementsDatabase[symbol]
        except KeyError:
            print 'Element: ', symbol, ' not present!'
            raise

        # most abundant mass number (if mass number is not provided)
        aux = AtomicElementsDatabase['isotopes_'+symbol]['most_abundant']
        self.data['mass'] = AtomicElementsDatabase['isotopes_'+symbol][aux]['atomicWeight']

        if not BOA:            
            if mass_number != 0:
                try:
                    self.data.update(AtomicElementsDatabase['isotopes_'+symbol][mass_number])
                    self.data['mass'] = self.data['atomicWeight']
                except KeyError:
                    print 'Mass number: ', str(mass_number), symbol, ' not present!'
                    raise
            else:
                self.data.update(AtomicElementsDatabase['isotopes_'+symbol][aux])

        #converting to Bohr
        position = np.array(position, dtype=np.float64)
        if units == 'Angstroms' : position *= UnitsDatabase['Bohr']
        self.data['position'] = position

    def isQuantum(self): 
        """ Returns whether the nuclei of the atomic element is being treated in the 
        BOA framework or not"""       
        try:
            self.data['mass_number']
            return True
        except KeyError:
            return False

    def getValue(self, key):
        """Returns the value stored in key
        """
        try:
            return self.data[key]
        except KeyError:
            raise
    
    def setValue(self, key, value):
        """Returns the value stored in key
        """
        self.data[key] = value

    def show(self):
        """Shows the information of the object
        """
        print self.data
