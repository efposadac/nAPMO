# file: test_stack.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it
from interfaces.stack import *

def test_stack_interface():
	a = Stack()
	assert a.size() == 0
	a.push(2)
	assert a.size() == 1
	assert a.isEmpty() == False
	assert a.pop() == 2 
	assert a.isEmpty() == True
	a.push(0)
	a.push(1)
	a.push(2)
	assert a.size() == 3
	for i in range(a.size()):
		assert a.get(i) == i
	try:
		a.get(100)
		assert False, 'Expecting Failure!'
	except:
		pass
	assert a.peek() == 2
	assert a.size() == 3

