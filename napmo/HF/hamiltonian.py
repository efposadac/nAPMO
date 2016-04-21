# file: hamiltonian.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co


class Hamiltonian(object):

    """
    docstring for Hamiltonian
    """

    def __init__(self, molecule):
        super(Hamiltonian, self).__init__()
        self.potential = np.zeros(molecule.basis)
        self.exchange = exchange
        self.kinetic = kinetic
        self.attraction

    def add_operator()
