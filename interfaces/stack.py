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
        """
        Clean the stack
        """
        return self.items == []

    def push(self, item):
        """Adds a new item in the stack
        """
        self.items.append(item)

    def pop(self):
        """Return the last element in the stack
        """
        return self.items.pop()

    def peek(self):
        """ Deletes the last element in the stack
        """
        return self.items[len(self.items)-1]

    def size(self):
        """ Returns the size of the stack
        """
        return len(self.items)

