# file: elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np


class Stack(object):
    """Defines a stack class as a list.
    """
    def __init__(self):
        super(Stack, self).__init__()

        self.items = []

    def isEmpty(self):
        """ Check whether the stack is empty or not.
        """
        return self.items == []

    def push(self, item):
        """Adds a new item in the stack
        """
        self.items.append(item)

    def pop(self):
        """Pop the last element in the stack
        """
        assert self.size() > 0
        return self.items.pop()

    def peek(self):
        """Return the last element in the stack
        """
        assert self.size() > 0
        return self.items[len(self.items)-1]

    def size(self):
        """ Returns the size of the stack
        """
        return len(self.items)

    def get(self, i):
        """ Returns the ith object in the stack
        """
        assert isinstance(i, int)
        assert i < self.size()
        return self.items[i]
