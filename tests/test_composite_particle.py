# file: test_composite_particle.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from interfaces.composite_particle import *

def test_composite_particle_interface():
	a = CompositeParticle()
	assert a.name() == 'none'
	
	a = CompositeParticle('H2O')
	assert a.name() == 'H2O'

	a.add_atom('C', [0.,0.,0.])
	assert a.size() == 2 # e- and C

	a.add_elementary_particle('e-',[0,0,0])
	assert a.size() == 2

	a.add_elementary_particle('u-',[0,0,0])

	assert a.size() == 3 #Keys in composite particle
	assert a.size('e-') == 2 # Size of stack e-
	assert a.n_particles('e-') == 7 # number of particles of kind e-

	try:
		a.n_particles('x-')
		assert False, 'Failure expected!'
	except:
		pass

	try:
		a.show()
		assert True
	except:
		pass

	try:
		a.get_particle('e+')
		assert False, 'Failure expected!'
	except:
		pass

	assert a.get_particle('u-').peek().get('name') == 'muon'

	assert a.get_particle('atoms', 0).isQuantum() == False

	a.add_atom('H', [0.,0.,0.], BOA=False)

	assert a.get_particle('atoms', 1).isQuantum() == True
