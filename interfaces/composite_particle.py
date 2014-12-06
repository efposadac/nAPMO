# file: CompositeParticle.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np
from copy import deepcopy

from interfaces.atomic_element import AtomicElement
from interfaces.elementary_particle import ElementaryParticle
from utilities.databases import UnitsDatabase
from interfaces.stack import Stack

class CompositeParticle(object):
	"""Defines a composite particle, i.e. atoms, molecules or nuclei.

		since a composite particle can be a molecule, you can add different atoms and 
		elementary particles to the molecule.
	"""
	def __init__(self):
		super(CompositeParticle, self).__init__()

		self.atoms = Stack()
		self.data = {}
		
	def add_atom(self, symbol, position, BOA=True, mass_number=0, units='Angstroms'):
		"""Adds an atom to the composite particle. 
		"""
		#converting to Bohr
		position = np.array(position, dtype=np.float64)
		if units == 'Angstroms' : position *= UnitsDatabase['Bohr']

		atom = {}
		atom[symbol] = AtomicElement(symbol, position=position, BOA=BOA, mass_number=mass_number, units='Bohr')
		self.atoms.push(deepcopy(atom))
		self.add_elementary_particle('e-', position, atom[symbol].get('atomicNumber'), units='Bohr')
		if not BOA: self.add_nuclei(atom[symbol])

	def add_nuclei(self, atom):
		"""Adds a 'quantum nuclei'
		"""
		symbol = str(atom.get('symbol'))+'_'+str(atom.get('mass_number'))

		try:
			self.data[symbol][len(self.data[symbol])-1].data['symbol']
		except KeyError:
			self.data[symbol] = []
		
		self.data[symbol].append(deepcopy(atom))
		self.data[symbol][len(self.data[symbol])-1].data['size'] = 1
		self.data[symbol][len(self.data[symbol])-1].data['symbol'] = symbol

	def add_elementary_particle(self, symbol, position, size=1, units='Angstroms'):
		"""Adds an elementary particle into the composite particle.
		"""
		#converting to Bohr
		position = np.array(position, dtype=np.float64)
		if units == 'Angstroms' : position *= UnitsDatabase['Bohr']

		try:
			self.data[symbol][len(self.data[symbol])-1].data['symbol']

		except KeyError:
			self.data[symbol] = []
		
		self.data[symbol].append(deepcopy(ElementaryParticle(symbol, position=position)))
		self.data[symbol][len(self.data[symbol])-1].data['size'] = size

	def show(self):
		"""Shows information about particles
		"""
		print ''
		print '** Quantum particles **'
		print ' '
		for key in self.data.iterkeys():			
			print 'Info:', key

		if self.atoms.size() > 0:
			print ''
			print '** Atoms **'
			print ' '
			print "Symbol     Position (Bohr)"
			print ''
			for i in xrange(self.atoms.size()):
				for item in self.atoms.get(i):
					atom = self.atoms.get(i)[item].get
					print '{0:10}'.format(atom('symbol')), atom('position')
