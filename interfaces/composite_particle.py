# file: CompositeParticle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from interfaces.atomic_element import AtomicElement
from interfaces.elementary_particle import ElementaryParticle
from utilities.databases import UnitsDatabase
from interfaces.stack import Stack


class CompositeParticle(object):
    """Defines a composite particle, i.e. atoms, molecules or nuclei.
        Since a composite particle can be a molecule, it can be added different atoms and
        elementary particles to the molecule.
    """
    def __init__(self, name='none'):
        super(CompositeParticle, self).__init__()

        self._name = name
        self.data = {}

    def add_atom(self, symbol, position, BOA=True, mass_number=0, units='Angstroms'):
        """Adds an atom to the composite particle.
        """
        assert isinstance(symbol, str)
        assert len(position) == 3
        assert isinstance(BOA, bool)
        assert isinstance(mass_number, int)
        assert isinstance(units, str)

        try:
            self.data['atoms']
        except KeyError:
            self.data['atoms'] = Stack()

        # Converting to Bohr
        position = np.array(position, dtype=np.float64)
        if units == 'Angstroms':
            position *= UnitsDatabase['Bohr']

        atom = AtomicElement(symbol, position=position, BOA=BOA, mass_number=mass_number, units='Bohr')
        self.data['atoms'].push(atom)
        self.add_elementary_particle('e-', position, size=atom.get('atomicNumber'), units='Bohr')

    def add_elementary_particle(self, symbol, position, size=1, units='Angstroms'):
        """Adds an elementary particle into the composite particle.
        """
        assert isinstance(symbol, str)
        assert len(position) == 3
        assert isinstance(size, int)
        assert isinstance(units, str)

        # Converting to Bohr
        position = np.array(position, dtype=np.float64)

        if units == 'Angstroms':
            position *= UnitsDatabase['Bohr']

        try:
            self.data[symbol].peek().set('symbol', symbol)

        except KeyError:
            self.data[symbol] = Stack()

        self.data[symbol].push(ElementaryParticle(symbol, position=position))
        self.data[symbol].peek().set('size', size)

    def name(self):
        """Returns the name.
        """
        return self._name

    def size(self, symbol='all'):
        """Returns the size of the composite particle.
        If symbol is given return the size of the stack 'symbol'.
        """

        assert isinstance(symbol, str)

        if symbol == 'all':
            aux = [i for i in self.data.iterkeys()]
            return len(aux)
        else:
            try:
                return self.data[symbol].size()
            except KeyError:
                raise

    def n_particles(self, symbol):
        """Returns the number of particles of kind 'symbol'.
        """

        assert isinstance(symbol, str)

        size = 0

        try:
            for i in xrange(self.size(symbol)):
                size += self.get_particle(symbol, i).get('size')
            return size
        except KeyError:
            raise

    def iter_symbol(self):
        """Return all symbols of (particles) in the composite object.
        """
        return self.data.iterkeys()

    def get_particle(self, symbol, iterator=-1):
        """Returns  particle 'symbol'
        """
        assert isinstance(symbol, str)

        try:
            if iterator == -1:
                return self.data[symbol]
            else:
                assert iterator < self.size(symbol)
                return self.get_particle(symbol).get(iterator)
        except KeyError:
            raise

    def show(self):
        """Shows information about particles
        """
        print '==============================================='
        print 'Object: '+type(self).__name__, self.name()
        for symbol in self.iter_symbol():
            if symbol != 'atoms':
                print type(self.get_particle(symbol).peek()).__name__, ': ', symbol,
                print '  n. particles: ', self.n_particles(symbol)
            else:
                print '-----------------------------------------------'
                print type(self.get_particle(symbol).peek()).__name__
                print "Symbol     Position (Bohr)"
                for i in xrange(self.size(symbol)):
                    atoms = self.get_particle(symbol).get(i).get
                    print '{0:10}'.format(atoms('symbol')), atoms('position')
                print '-----------------------------------------------'
