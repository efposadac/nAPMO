# file: elementary_particle.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
import numpy as np


class Stack(list):
    """
    Defines a stack class as a list.
    """
    def __init__(self):
        super(Stack, self).__init__()

    def isEmpty(self):
        """
        Check whether the stack is empty or not.

        Returns:
            bool: Whether the stack is empty or not.
        """
        return not self

    def push(self, item):
        """
        Adds a new item in the stack

        Args:
            item (any): Item to be added in the stack.
        """
        self.append(item)

    def peek(self):
        """
        Returns:
            any: The last element in the stack
        """
        assert len(self) > 0
        return self[-1]
