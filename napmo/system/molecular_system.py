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

from napmo.system.atomic_element import AtomicElement
from napmo.system.elementary_particle import ElementaryParticle
from napmo.system.basis_set import BasisSet, BasisSet_C
from napmo.data.constants import ANGSTROM_TO_BOHR


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
            origin (ndarray): Origin of the atom (Cartesian coordinates)
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

        atom = AtomicElement(symbol, origin=origin, BOA=BOA,
                             mass_number=mass_number, units='Bohr')

        if basis_name is None:
            basis_name = basis_file

        self.add_elementary_particle(
            'e-', size=atom.get('atomic_number'), units='BOHR', basis_name=basis_name)

        # load basis-set
        atom['basis'] = BasisSet(basis_name)

        if basis_file is not None:
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                atom.get('basis').load_gaussian(
                    symbol, basis_data, origin)
                self.get('e-')['basis'].load_gaussian(symbol,
                                                      basis_data, origin)

            elif basis_kind == 'STO':
                atom.get('basis').load_slater(
                    symbol, basis_data, origin)
                self.get('e-')['basis'].load_slater(
                    symbol, basis_data, origin)

        self.get('atoms').append(atom)

    def add_elementary_particle(self, symbol,
                                origin=None, size=1, units='ANGSTROMS',
                                basis_kind='GTO', basis_name=None, basis_file=None):
        """
        Adds an elementary particle into the molecular system.

        Args:
            symbol (str): Symbol of the elementary particle.
            origin (ndarray): Origin of the elementary particle (Cartesian coordinates)
            size (int): Number of elementary particle to be added.
            units (str, optional): Units of the origin, valid values are 'ANGSTROMS' or 'Bohr'
        """
        assert isinstance(symbol, str)
        assert isinstance(size, int)
        assert isinstance(units, str)

        if symbol not in self:
            self[symbol] = ElementaryParticle(symbol)
            self[symbol]['size'] = 0
            self[symbol]['basis'] = BasisSet(basis_name)
            self[symbol]['origin'] = []

        # Converting to Bohr
        if origin is not None:
            assert len(origin) == 3
            origin = np.array(origin, dtype=np.float64)

            if units == 'ANGSTROMS':
                origin *= ANGSTROM_TO_BOHR

            self[symbol]['origin'].append(origin)

        self[symbol]['size'] += size

        # load basis-set
        if basis_file is not None:
            file = open(basis_file)
            basis_data = file.read().replace('\n', '')
            file.close()

            if basis_kind == 'GTO':
                self[symbol].get('basis').load_gaussian(
                    symbol, basis_data, origin)
            elif basis_kind == 'STO':
                self[symbol].get('basis').load_slater(
                    symbol, basis_data, origin)

    def n_occupation(self, symbol):
        """
        Returns the occupation number for particle ``symbol``.

        Returns:
            int: Occupation of quantum species ``symbol``.
        """
        return ceil(self.n_particles(symbol) * self[symbol]['spin'])

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
            if symbol != 'atoms':
                output = self[symbol].get('size')

        return output

    def get_basis(self, symbol):
        """
        Returns the basis set of the system of a given ``symbol`` particle.

        Returns:
            BasisSet: basis set
        """
        basis = None

        if symbol in self:
            basis = self[symbol].get('basis')
        else:
            for atom in self.get('atoms'):
                if atom.get('symbol') == symbol:
                    basis = atom.get('basis')

        if isinstance(basis, BasisSet):
            return basis
        else:
            raise KeyError('The symbol:', symbol, 'does not exist in the molecular system!')

    def get_basis_as_cstruct(self, symbol):
        """
        Returns the basis set of the system of a given ``symbol`` particle as a CTYPES struct.

        Returns:
            BasisSet_C: basis set CTYPES struct
        """
        basis = self.get_basis(symbol)
        return BasisSet_C(basis)

    def show(self):
        """
        Shows information about particles
        """
        print('==================================================================')
        print('Object: ' + type(self).__name__)
        print('------------------------------------------------------------------')
        for symbol in self.keys():
            particles = self.get(symbol)
            if symbol is 'atoms':
                print(type(particles[-1]).__name__, ': ', symbol, end=" ")
                print(' number of e-: ', self.n_particles('e-'))
                print(
                    '{0:7}'.format("Symbol"),
                    '{0:5}'.format("Z"),
                    '{0:40}'.format("origin (Bohr)"),
                    '{0:10}'.format("Basis-set")
                )
                for particle in particles:
                    print(
                        '{0:7}'.format(str(particle.get('symbol'))),
                        '{0:5}'.format(str(particle.get('atomic_number'))),
                        '{0:40}'.format(str(particle.get('origin'))),
                        '{0:10}'.format(str(particle.get('basis')['name']))
                    )
            else:
                if len(particles.get('origin')) > 1:
                    print(type(particles).__name__, ': ', symbol, end=" ")
                    print(' number of particles: ', self.n_particles(symbol))
                    print(
                        '{0:7}'.format("Symbol"),
                        '{0:5}'.format("N"),
                        '{0:40}'.format("origin (Bohr)"),
                        '{0:10}'.format("Basis-set")
                    )
                for i in range(len(particles.get('origin'))):
                    print(
                        '{0:7}'.format(str(particles.get('symbol'))),
                        '{0:5}'.format(str(particles.get('size'))),
                        '{0:40}'.format(str(particles.get('origin')[i])),
                        '{0:10}'.format(str(particles.get('basis')['name']))
                    )
        print('------------------------------------------------------------------')
