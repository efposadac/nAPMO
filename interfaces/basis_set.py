# file: basis_set.py
# nAPMO package
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
import numpy as np
import json

from interfaces.stack import Stack
from interfaces.contracted_slater import ContractedSlater
from interfaces.contracted_gaussian import ContractedGaussian


class BasisSet(dict):
    """
    Basis-set interface. (dict)

    This class allows the management of STO and GTO basis-sets.
    """
    def __init__(self, name='user'):
        super(BasisSet, self).__init__()
        self['name'] = name
        self['particle'] = None
        self['function'] = Stack()
        self['kind'] = None
        self['json'] = None
        self['length'] = 0

    def __add__(self, other):
        self['name'] += ', '+other.get('name')
        self['kind'] += ', '+other.get('kind')
        self['particle'] += ', '+other.get('particle')
        self['function'] += other.get('function')
        self['json'] += other.get('json')
        self['length'] += other.get('length')
        self['t_length'] += other.get('t_length')

        return self

    def load_gaussian(self, particle, data, origin=[0.0, 0.0, 0.0]):
        """
        Load a Gaussian Type Orbital (GTO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self['particle'] = particle
        self['kind'] = 'GTO'
        self['json'] = json.loads(data)[particle]

        lvalue = {'s': 0, 'p': 1, 'd': 2, 'f': 3, 'g': 4}

        for k in range(len(self['json'])):
            l = lvalue[self['json'][k]['angular']]
            for i in range(l+1):
                x = l - i
                for j in range(i+1):
                    y = i - j
                    z = j
                    self['function'].push(ContractedGaussian(
                        np.array(self['json'][k]['prim']),
                        np.array(self['json'][k]['cont']),
                        np.array(origin),
                        np.array([x, y, z])
                        ))

        self['length'] = len(self.get('function'))

        self['t_length'] = 0
        for function in self.get('function'):
            self['t_length'] += function.get('length')

    def load_slater(self, particle, data, origin=[0.0, 0.0, 0.0]):
        """
        Load a Slater Type Orbital (STO) basis-set.

        Args:
            particle (str): Symbol of the particle.
            data (json): json formatted data to be loaded.
        """
        self['particle'] = particle
        self["kind"] = 'STO'
        self["json"] = json.loads(data)[particle]

        lvalue = {"s": 0, "p": 1, "d": 2, "f": 3, "g": 4}

        for k in range(len(self["json"])):
            l = lvalue[self["json"][k]["angular"]]
            for m in range(-l, l+1):
                self["function"].append(ContractedSlater(
                    np.array(self["json"][k]["prim"]),
                    np.array(self["json"][k]["cont"]),
                    np.array(origin),
                    np.array(self["json"][k]["n"]),
                    l,
                    m
                ))
        self['length'] = len(self.get('function'))

        self['t_length'] = 0
        for function in self.get('function'):
            self['t_length'] += function.get('length')

    def compute(self, coord=np.array([0.0, 0.0, 0.0])):
        """
        Compute all the basis-set functions at given ``coord``.

        Args:
            coord (numpy.ndarray(3)): coordinates in which the basis set will be evaluated.

        """
        output = np.zeros(self.get('length'))
        for i in range(self.get('length')):
            output[i] = self.get('function')[i].compute(coord)

        return output

    def show_json(self):
        """
        Prints the basis-set in json format.
        """
        print(json.dumps(self["json"], sort_keys=True, indent=4))

    def show(self):
        """
        Prints extended information of the basis-set object.
        """
        print("================")
        print("Basis set info")
        print("================")
        print("Name: ", self.get('name'))
        print("Particle: ", self.get('particle'))
        print("Kind: ", self.get('kind'))
        print("Length: ", self.get('length'))
        print("****************")
        print("Functions info: ")
        print("****************")
        i = 1
        for function in self.get('function'):
            print("")
            print("*** Function: ", i)
            print("")
            function.show()
            i += 1
