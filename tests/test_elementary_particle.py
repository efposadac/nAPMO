# file: test_elementary_particle.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from interfaces.elementary_particle import *

def test_elementary_particle_interface():
	a = ElementaryParticle()
	assert a.get('name') == 'user'
	a.set('name', 'electron')
	assert a.get('name') == 'electron'
	a = ElementaryParticle('e-',position=[0,1,2])
	assert a.get('name') == 'electron'
	assert a.get('symbol') == 'e-'
	assert a.get('charge') == -1.
	assert a.get('mass') == 1.
	for i in range(3):
		assert a.get('position')[i] == float(i)
	try:
		a.get('address')
		assert False, 'Failure expected!'
	except:
		pass
	try:
		a.show()
		assert True
	except:
		pass