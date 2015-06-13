# file: molecular_system.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from interfaces.atomic_element import AtomicElement
from interfaces.elementary_particle import ElementaryParticle
from utilities.constants import ANGSTROM_TO_BOHR
from interfaces.stack import Stack
from interfaces.basis_set import BasisSet


class MolecularSystem(dict):
    """
    Defines a molecular system, containing different kinds of quantum 'species' i.e. atoms, muons, positrons, etc.

    Args:
        name (str): Name of the object.
    """
    def __init__(self):
        super(MolecularSystem, self).__init__()

    def add_atom(self, symbol, origin, BOA=True, mass_number=0, units='Angstroms',
                 basis_kind='GTO', basis_name='none', basis_file='none'):
        """
        Adds an atom to the molecular system.

        Args:
            symbol (str): Symbol of the atom.
            origin (array[3]): Origin of the atom (Cartesian coordinates)
            BOA (bool, optional): Whether the atom nuclei will be treated in the BOA approach or not. Default is True
            mass_number (int, optional): Mass number of element 'symbol' :math:`:= A = p^+ + n^o`.
                If 0 the system will choose the most abundant one. Default is 0.
            units (str, optional): Units of the origin, valid values are 'Angstroms' or 'Bohr'
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(BOA, bool)
        assert isinstance(mass_number, int)
        assert isinstance(units, str)

        if 'atoms' not in self:
            self['atoms'] = Stack()

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'Angstroms':
            origin *= ANGSTROM_TO_BOHR

        atom = AtomicElement(symbol, origin=origin, BOA=BOA, mass_number=mass_number, units='Bohr')
        atom['size'] = 1

        # load basis-set
        atom['basis'] = BasisSet(symbol, origin, basis_name)

        if basis_file != 'none':
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                atom.get('basis').load_gaussian(basis_data)
            elif basis_kind == 'STO':
                atom.get('basis').load_slater(basis_data)

        self['atoms'].push(atom)
        self.add_elementary_particle('e-', origin, size=atom.get('atomic_number'), units='Bohr')

    def add_elementary_particle(self, symbol, origin, size=1, units='Angstroms',
                                basis_kind='GTO', basis_name='none', basis_file='none'):
        """
        Adds an elementary particle into the molecular system.

        Args:
            symbol (str): Symbol of the elementary particle.
            origin (array[3]): Origin of the elementary particle (Cartesian coordinates)
            size (int): Number of elementary particle to be added.
            units (str, optional): Units of the origin, valid values are 'Angstroms' or 'Bohr'
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(size, int)
        assert isinstance(units, str)

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)

        if units == 'Angstroms':
            origin *= ANGSTROM_TO_BOHR

        if symbol not in self:
            self[symbol] = Stack()

        particle = ElementaryParticle(symbol, origin=origin)

        self[symbol].push(particle)
        self[symbol].peek()['size'] = size

        # load basis-set
        self[symbol].peek()['basis'] = BasisSet(symbol, origin, basis_name)

        if basis_file != 'none':
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                self[symbol].peek().get('basis').load_gaussian(basis_data)
            elif basis_kind == 'STO':
                self[symbol].peek().get('basis').load_slater(basis_data)

    def n_elementary_particles(self):
        """
        Calculates the number of elementary particles in the system.

        Returns:
            int: Number of quantum species in the system.
        """
        output = len(self)
        if 'atoms' in self:
            output -= 1

        return output

    def n_atoms(self):
        """
        Calculates the number of atoms in the system

        Returns:
            int: Number of atoms in the object.
        """
        output = 0
        if 'atoms' in self:
            output = len(self['atoms'])

        return output

    def n_particles(self, symbol):
        """
        Return the number of particles of a given ``symbol`` if this is present in the object.

        Returns:
            int: Number of particles of symbol ``symbol``
        """
        output = 0
        if symbol in self:
            for item in self[symbol]:
                output += item.get('size')

        return output

    def show(self):
        """
        Shows information about particles
        """
        print('===============================================')
        print('Object: ' + type(self).__name__)
        for symbol in self.keys():
            if symbol != 'atoms':
                particle = self.get(symbol)
                print(type(particle.peek()).__name__, ': ', symbol, end=" ")
                print('  n. particles: ', self.n_particles(symbol))
            else:
                particle = self.get('atoms')
                print('-----------------------------------------------')
                print(type(particle.peek()).__name__)
                print("Symbol     origin (Bohr)                       Basis-set")
                for i in particle:
                    atoms = i.get
                    print('{0:10}'.format(atoms('symbol')), atoms('origin'), atoms('basis')['name'])
                print('-----------------------------------------------')
