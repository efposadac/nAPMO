# file: molecular_system.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np
import copy

import napmo


class MolecularSystem(dict):

    """
    Defines a molecular system, containing different kinds of quantum 'species'
    i.e. atoms, muon, positrons, etc. (dict)

    Args:
        name (str): Name of the object.
    """

    def __init__(self):
        super(MolecularSystem, self).__init__()
        self._point_charges = []

    def add_atom(self, symbol, origin, units='ANGSTROMS', quantum=False,
                 basis_name=None, basis_file=None):
        """
        Adds an atom to the molecular system.

        Args:
            symbol (str): Symbol of the atom.
            origin (ndarray): Origin of the atom (Cartesian coordinates)
            units (str, optional): Units of the origin, valid values are
                'ANGSTROMS' or 'BOHR'
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(units, str)

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units is 'ANGSTROMS':
            origin *= napmo.ANGSTROM_TO_BOHR

        # Fetching atom database
        atom = napmo.AtomicElement(symbol, origin, 'BOHR')

        # load basis-set
        basis = None
        if basis_name:
            basis = napmo.BasisSet(
                basis_name, symbol, origin=origin, basis_file=basis_file)

            atom['basis'] = basis

        atom['is_quantum'] = quantum

        self.add_elementary_particle('e-', origin, units='BOHR', size=atom.get('atomic_number'),
                                     basis_name=basis_name, basis_file=basis_file,
                                     basis_set=basis, particle=atom)

    def add_elementary_particle(self, symbol, origin, units='ANGSTROMS', size=1,
                                basis_name=None, basis_file=None, basis_set=None,
                                particle=None):
        """
        Adds an elementary particle into the molecular system.

        Args:
            symbol (str): Symbol of the elementary particle.
            origin (ndarray): Origin of the elementary particle
                (Cartesian coordinates)
            units (str, optional): Units of the origin, valid values are
                'ANGSTROMS' or 'Bohr'
            size (int): Number of elementary particle to be added.
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(size, int)
        assert isinstance(units, str)

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units is 'ANGSTROMS':
            origin *= napmo.ANGSTROM_TO_BOHR

        # Fetching information from data
        eparticle = napmo.ElementaryParticle(symbol, origin, 'BOHR')

        # setting defaults for species
        self.setdefault(symbol, eparticle)

        self[symbol].setdefault('id', self.size_species - 1)
        self[symbol].setdefault('size', 0)
        self[symbol].setdefault('particles', [])
        self[symbol].setdefault('origin', [])
        self[symbol].pop('origin')

        # load basis-set
        basis = None
        if basis_name and not basis_set:
            basis = napmo.BasisSet(basis_name, symbol, origin=origin,
                                   basis_file=basis_file)
        if basis_set:
            basis = basis_set

        if 'basis' not in self[symbol] and basis:
            self.get(symbol)['basis'] = napmo.BasisSet(basis_name, symbol)

        if basis:
            self.get(symbol)['basis'].update(basis)
        else:
            eparticle['is_quantum'] = False

        self.get(symbol)['size'] += size
        self.get(symbol)['occupation'] = int(
            self.get(symbol)['size'] * self.get(symbol)['particlesfraction'])

        if particle:
            self.get(symbol)['particles'].append(particle)
        else:
            self.get(symbol)['particles'].append(eparticle)

        if not self.get(symbol).get('particles')[-1].is_quantum:
            self._point_charges.append(self.get(symbol).get('particles')[-1])

    def add_nucleus(self, symbol, origin, units='ANGSTROMS', size=1,
                    basis_name=None, basis_file=None):
        """
        Adds a nucleus to the molecular system.

        Args:
            symbol (str): Symbol of the nucleus ie. H_2 (Deuterium).
            origin (ndarray): Origin of the nucleus (Cartesian coordinates)
            units (str, optional): Units of the origin, valid values are
                'ANGSTROMS' or 'BOHR'
            size (int): Number of particles to be added.
        """
        assert isinstance(symbol, str)
        assert len(origin) == 3
        assert isinstance(units, str)

        # Being sure this is a nucleus specification, otherwise, add an atom
        if symbol.find("_") < 0:
            self.add_atom(symbol, origin, units=units, quantum=False,
                          basis_name=basis_name, basis_file=basis_file)
            return

        # Converting to Bohr
        origin = np.array(origin, dtype=np.float64)
        if units is 'ANGSTROMS':
            origin *= napmo.ANGSTROM_TO_BOHR

        # Fetching atom object
        nucleus = napmo.AtomicElement(symbol, origin, 'BOHR')

        # setting defaults for species
        self.setdefault(symbol, {})

        self[symbol].setdefault('id', self.size_species - 1)
        self[symbol].setdefault('size', 0)
        self[symbol].setdefault('particles', [])
        self[symbol].setdefault('origin', [])
        self[symbol].pop('origin')

        # load basis-set
        basis = None
        if basis_name:
            basis = napmo.BasisSet(basis_name, symbol, origin, basis_file)

        if 'basis' not in self[symbol] and basis:
            self.get(symbol)['basis'] = napmo.BasisSet(basis_name, symbol)

        if basis:
            self.get(symbol)['basis'].update(basis)
            nucleus['basis'] = basis
        else:
            nucleus['is_quantum'] = False

        # Adding coupling constants
        self.get(symbol).update(
            napmo.CouplingConstantsDatabase()[symbol.lower()])

        # Add to the molecular system
        self.get(symbol)['size'] += size
        self.get(symbol)['occupation'] = int(
            self.get(symbol)['size'] * self.get(symbol)['particlesfraction'])
        self.get(symbol)['charge'] = nucleus.get('atomic_number', 1)
        self.get(symbol)['spin'] = nucleus.get('spin', 1)
        self.get(symbol)['mass'] = nucleus.get('mass', 1)
        self.get(symbol)['particles'].append(nucleus)

        if not self.get(symbol).get('particles')[-1].is_quantum:
            self._point_charges.append(self.get(symbol).get('particles')[-1])

    def set_charges(self, data, open_shell=False):

        assert isinstance(data, dict)

        for key in data:
            if key not in self:
                print("Imposible to set " + key +
                      " charges: Particle does not exist!")
                return

            multi = data.get(key, {}).get('multiplicity', None)
            charge = data.get(key, {}).get('charge', 0)

            np = self.size_particles(key) + charge

            # electronic case
            if key == 'e-':
                if not multi:
                    multi = 2 * ((np % 2) * 0.5) + 1

                spin = (multi - 1) * 0.5

                nereq = spin * 2
                eleft = np - nereq

                if eleft % 2 != 0:
                    print("Bad charge / multiplicity")
                    raise ValueError

                alpha = (np - np % 2) * 0.5
                alpha += spin / napmo.SPIN_ELECTRON
                beta = np - alpha
                keys = {'e-alpha': alpha, 'e-beta': beta}

                if alpha != beta or open_shell:
                    for k in keys:
                        eparticle = napmo.ElementaryParticle(k)
                        eparticle.pop('origin')

                        self[k] = copy.deepcopy(self.get('e-'))

                        self.get(k).update(eparticle)
                        self.get(k)['size'] = keys.get(k)

                        self.get(k)['occupation'] = int(
                            self.get(k)['size'] * self.get(k)['particlesfraction'])

                    self.pop('e-')

                    for i, k in enumerate(self):
                        self.get(k)['id'] = i

                else:
                    self.get('e-')['size'] += charge
                    self.get('e-')['occupation'] = int(
                        self.get('e-')['size'] * self.get('e-')['particlesfraction'])
            else:
                self.get(key, {})['size'] += charge
                self.get(key, {})['occupation'] = int(
                    self.get(key, {})['size'] * self.get(key, {})['particlesfraction'])

    def size_particles(self, symbol):
        """
        Returns the number of particles of a given species ``symbol`` in the object.

        Returns:
            int: Number of particles of the species ``symbol``
        """
        return self.get(symbol, {}).get('size', 0)

    def get_basis(self, symbol):
        """
        Returns the basis set of the system of a given ``symbol`` particle.

        Returns:
            napmo.BasisSet: basis set
        """
        out = self.get(symbol, {}).get('basis', None)

        if not out:
            for species in self:
                for particle in self.get(species, {}).get('particles', []):
                    if particle.get('symbol') is symbol:
                        out = particle.get('basis', None)
        return out

    def get_species(self, sid):
        """
        Return a list with the species ``symbol`` or `sid`

        Args:
            symbol(str): symbol of the species. ie. ``e-``, ``e+``, etc. Default ``None``
            sid(int): the species if. sid > 0. Default ``None``
        """
        for species in self:
            if self.get(species, {}).get('id', -100) is sid:
                return self.get(species)

    def get_basis_as_cstruct(self, symbol):
        """
        Returns the basis set of the system of a given ``symbol`` particle as a
        CTYPES structure.

        Returns:
            napmo.BasisSet_C: basis set CTYPES structure
        """
        out = None
        basis = self.get_basis(symbol)
        if basis:
            out = napmo.BasisSet_C(basis)
        return out

    @property
    def size_species(self):
        """
        Returns the number of quantum species in the system.

        Returns:
            int: Number of quantum species in the system.
        """
        return len(self)

    @property
    def point_charges(self):
        return self._point_charges

    def __repr__(self):

        out = """
==================================================
Object: {0:9s}
--------------------------------------------------

{1:7s} {2:4s} {3:^29s} {4:9s}
--------------------------------------------------
""".format(
            type(self).__name__,
            "Symbol",
            "N",
            "origin(BOHR)",
            "Basis"
        )

        basis = set()

        for species in self.keys():
            if species == 'e-beta':
                continue
            particles = self.get(species)
            for particle in particles.get('particles', []):

                out += '{0:3s} {1:5d} {2:10.6f} {3:10.6f} {4:10.6f} {5:10s} \n'.format(
                    particle.get('symbol', '--'),
                    particle.get('atomic_number', particle.get('size', 0)),
                    particle.get('origin', [0.0, 0.0, 0.0])[0],
                    particle.get('origin', [0.0, 0.0, 0.0])[1],
                    particle.get('origin', [0.0, 0.0, 0.0])[2],
                    particle.get('basis', {}).get('name', '--'))

                basis.add(str(particle.get('basis', '--')))

        out += ('--------------------------------------------------')

        return out + ''.join(s for s in basis)
