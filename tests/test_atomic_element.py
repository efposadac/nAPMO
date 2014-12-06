# file: test_atomic_element.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from interfaces.atomic_element import AtomicElement
import numpy as np

def test_atomic_element_interface():
	try:
		a = AtomicElement()
		assert False, 'Failure expected!'
	except:
		pass
	try:
		a = AtomicElement('XXX')
		assert False, 'Failure expected'
	except:
		pass
	a = AtomicElement('C')
	assert a.get('name') == 'Carbon'
	assert a.isQuantum() == False
	for i in xrange(3):
		assert a.get('position')[i] == 0.

	a = AtomicElement('H', BOA=False)

	assert a.isQuantum() == True
	assert a.get('mass_number') == 1

	a.set('mass_number',2)
	assert a.get('mass_number') == 1

	a = AtomicElement('H', mass_number = 3)
	assert a.isQuantum() == True

	try:
		a = AtomicElement('H', mass_number = 100)
		assert False, 'Failure expected'
	except:
		pass

	try:
		a.get('address')
		assert False, 'Failure expected'
	except :
		pass

	a.set('address','sky')
	assert a.get('address') == 'sky'
	
	try:
		a.show()
		assert True
	except:
		pass

