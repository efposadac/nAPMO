# file: elementary_particle.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from databases import ElementaryParticlesDatabase

class ElementaryParticle(object):
	"""An abstract python interface to create an elementary quantum particle
	i.e Leptons as electron, muon, etc

	Database information from:
	From: http://physics.nist.gov/constants
	"""
	def __init__(self, symbol='null', position=[0.0,0.0,0.0]):

		"""Pure interface class. It generates an "generic" particle by default.
		"""
		super(ElementaryParticle, self).__init__()

		try:
			self.data = ElementaryParticlesDatabase[symbol]
		except KeyError:
			print 'Elementary particle: ', symbol, ' not present!, creating one.'
			self.data = ElementaryParticlesDatabase['user'] 

		self.data['position'] = np.array(position, dtype=np.float64)

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
