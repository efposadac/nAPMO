# file: test_stack.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import os
import sys

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.stack import *


def test_stack_interface():
    a = Stack()
    assert len(a) == 0
    a.push(2)
    assert len(a) == 1
    assert a.isEmpty() is False
    assert a.pop() == 2
    assert a.isEmpty() is True
    a.push(0)
    a.push(1)
    a.push(2)
    assert len(a) == 3

    for i in range(len(a)):
        assert a[i] == i

    try:
        a[100]
        assert False, 'Expecting Failure!'
    except:
        assert True

    assert a.peek() == 2
    assert len(a) == 3

# test_stack_interface()
