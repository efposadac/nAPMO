# file: molecular_system.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

from copy import deepcopy

import numpy as np
from math import ceil

from napmo.interfaces.atomic_element import AtomicElement
from napmo.interfaces.elementary_particle import ElementaryParticle
from napmo.utilities.constants import ANGSTROM_TO_BOHR
from napmo.interfaces.basis_set import BasisSet


class MolecularSystem(dict):
    """
    Defines a molecular system, containing different kinds of quantum 'species' i.e. atoms, muons, positrons, etc. (dict)

    Args:
        name (str): Name of the object.
    """
    def __init__(self):
        super(MolecularSystem, self).__init__()

    def add_atom(self, symbol, origin, BOA=True, mass_number=0, units='ANGSTROMS',
                 basis_kind='GTO', basis_name=None, basis_file=None):
        """
        Adds an atom to the molecular system.

        Args:
            symbol (str): Symbol of the atom.
            origin (numpy.ndarray(3)): Origin of the atom (Cartesian coordinates)
            BOA (bool, optional): Whether the atom nuclei will be treated in the BOA approach or not. Default is True
            mass_number (int, optional): Mass number of element 'symbol' :math:`:= A = p^+ + n^o`.
                If 0 the system will choose the most abundant one. Default is 0.
            units (str, optional): Units of the origin, valid values are 'ANGSTROMS' or 'BOHR'
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(BOA, bool)
        assert isinstance(mass_number, int)
        assert isinstance(units, str)

        if 'atoms' not in self:
            self['atoms'] = []

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units == 'ANGSTROMS':
            origin *= ANGSTROM_TO_BOHR

        atom = AtomicElement(symbol, origin=origin, BOA=BOA, mass_number=mass_number, units='Bohr')
        atom['size'] = 1

        # load basis-set
        if basis_name is None:
            basis_name = basis_file

        atom['basis'] = BasisSet(basis_name)

        if basis_file is not None:
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                atom.get('basis').load_gaussian(symbol, basis_data, origin)
            elif basis_kind == 'STO':
                atom.get('basis').load_slater(symbol, basis_data, origin)

        self.get('atoms').append(atom)
        self.add_elementary_particle('e-', origin, size=atom.get('atomic_number'), units='Bohr')
        self.get('e-')[-1]['basis'] = self.get('atoms')[-1].get('basis')

    def add_elementary_particle(self, symbol, origin, size=1, units='ANGSTROMS',
                                basis_kind='GTO', basis_name=None, basis_file=None):
        """
        Adds an elementary particle into the molecular system.

        Args:
            symbol (str): Symbol of the elementary particle.
            origin (array[3]): Origin of the elementary particle (Cartesian coordinates)
            size (int): Number of elementary particle to be added.
            units (str, optional): Units of the origin, valid values are 'ANGSTROMS' or 'Bohr'
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(size, int)
        assert isinstance(units, str)

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)

        if units == 'ANGSTROMS':
            origin *= ANGSTROM_TO_BOHR

        if symbol not in self:
            self[symbol] = []

        particle = ElementaryParticle(symbol, origin=origin)

        self[symbol].append(particle)
        self[symbol][-1]['size'] = size

        # load basis-set
        self[symbol][-1]['basis'] = BasisSet(basis_name)

        if basis_file != None:
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                self[symbol][-1].get('basis').load_gaussian(symbol, basis_data, origin)
            elif basis_kind == 'STO':
                self[symbol][-1].get('basis').load_slater(symbol, basis_data, origin)

    def n_occupation(self, symbol):
        """
        Returns the occupation number for particle ``symbol``.

        Returns:
            int: Occupation of quantum species ``symbol``.
        """
        return ceil(self.n_particles(symbol) * self[symbol][-1]['spin'])

    def n_elementary_particles(self):
        """
        Returns the number of elementary particles in the system.

        Returns:
            int: Number of quantum species in the system.
        """
        output = len(self)
        if 'atoms' in self:
            output -= 1

        return output

    def n_atoms(self):
        """
        Returns the number of atoms in the system

        Returns:
            int: Number of atoms in the object.
        """
        output = 0
        if 'atoms' in self:
            output = len(self['atoms'])

        return output

    def n_particles(self, symbol):
        """
        Returns the number of particles of a given ``symbol`` in the object.

        Returns:
            int: Number of particles of symbol ``symbol``
        """
        output = 0
        if symbol in self:
            for item in self[symbol]:
                output += item.get('size')

        return output

    def get_basis_set(self, symbol):
        if symbol in self:
            basis = deepcopy(self[symbol][0].get('basis'))
            for i in range(1, len(self[symbol])):
                basis += self[symbol][i].get('basis')
        else:
            for atom in self.get('atoms'):
                if atom.get('symbol') == symbol:
                    basis = atom.get('basis')

        if isinstance(basis, BasisSet):
            return basis
        else:
            raise KeyError

    def show(self):
        """
        Shows information about particles
        """
        print('==================================================================')
        print('Object: ' + type(self).__name__)
        print('------------------------------------------------------------------')
        for symbol in self.keys():
            particles = self.get(symbol)
            print(type(particles[-1]).__name__, ': ', symbol, end=" ")
            print('  n. particles: ', self.n_particles(symbol))
            print(
                    '{0:7}'.format("Symbol"),
                    '{0:5}'.format("n."),
                    '{0:40}'.format("origin (Bohr)"),
                    '{0:10}'.format("Basis-set")
                )
            for particle in particles:
                print(
                        '{0:7}'.format(particle.get('symbol')),
                        '{0:5}'.format(str(particle.get('size'))),
                        '{0:40}'.format(str(particle.get('origin'))),
                        '{0:10}'.format(str(particle.get('basis')['name']))
                    )
            print('------------------------------------------------------------------')
