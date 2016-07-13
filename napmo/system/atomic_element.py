# file: atomic_element.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np

import napmo


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
            'ANGSTROMS' or 'BOHR'
    '''

    def __init__(self, symbol, origin=np.zeros(3, dtype=np.float64), units='ANGSTROMS'):
        super(AtomicElement, self).__init__()

        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(units, str)

        data = napmo.AtomicElementsDatabase()

        _symbol = symbol
        mass_number = ''

        nuclei = symbol.find("_")
        if nuclei > 0:
            _symbol = symbol[:nuclei]
            mass_number = symbol[::-1][:nuclei]

        try:
            self.update(data[_symbol])
        except KeyError:
            print('Element: ', _symbol, ' not present!')
            raise

        # choose the most abundant mass number (if mass number is not provided)
        aux = data['isotopes_' + _symbol]['most_abundant']
        self['mass'] = (float(aux) * napmo.NEUTRON_MASS +
                        self['atomic_number'] *
                        (napmo.PROTON_MASS - napmo.NEUTRON_MASS))

        if mass_number:
            try:
                self.update(data['isotopes_' + _symbol][mass_number])
                self['mass'] = (int(mass_number) * napmo.NEUTRON_MASS +
                                self['atomic_number'] *
                                (napmo.PROTON_MASS - napmo.NEUTRON_MASS))

            except KeyError:
                print('Isotope: ', symbol, ' not present!')
                raise

        # converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'ANGSTROMS':
            origin *= napmo.ANGSTROM_TO_BOHR

        self['origin'] = origin
        self['symbol'] = symbol

    @property
    def is_quantum(self):
        """
        Returns:
            bool: Whether the nuclei of the atomic element is being treated in
                the BOA framework
        """
        return self.get('is_quantum')

    def __repr__(self):

        out = """
==================================================
Object:  {0:9s}
--------------------------------------------------
Name:    {1:9s}
Symbol:  {2:9s}
Quantum: {3:9s}
Z:       {4:<4d}
Mass:    {5:<5.3f}
Origin:  {6:<5.3f} {7:<5.3f} {8:<5.3f}
Basis:   {9:9s}
--------------------------------------------------""".format(
            type(self).__name__,
            self.get('name'),
            self.get('symbol'),
            str(self.is_quantum),
            self.get('atomic_number'),
            self.get('mass'),
            self.get('origin')[0],
            self.get('origin')[1],
            self.get('origin')[2],
            self.get('basis', {}).get('name', "None"),
        )

        out += ''.join(str(self.get('basis', {})))
        return out
