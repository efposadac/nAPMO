# file: napmo.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import print_function

import napmo
import matplotlib.pylab as plt

import sys


class NAPMO(object):
    """
    NAPMO system manager

    Args:
        data (InputParser) : parsed input data
    """

    Energy = 0.0

    def __init__(self, data, pprint=True):
        super(NAPMO, self).__init__()

        assert isinstance(data, napmo.InputParser)

        self._energy = 0.0

        self.data = data

        # Build system
        self.build_system()

        # Show system info
        if pprint:
            self.show()

    def build_system(self):
        """
        Builds the MolecularSystem object for the calculation
        """
        self.system = napmo.MolecularSystem()

        atomic_data = napmo.AtomicElementsDatabase()
        system_data = self.data.data

        for symbol in system_data:
            for p in system_data[symbol]:

                if symbol in atomic_data:  # load atoms
                    self.system.add_atom(
                        symbol, p.get('origin'),
                        basis_name=p.get('basis_name'),
                        basis_file=p.get('basis_file'),
                        quantum=p.get('quantum', False))

                elif symbol.find("_") > 0:  # load nuclei
                    self.system.add_nucleus(
                        symbol, p.get('origin'),
                        basis_name=p.get('basis_name'),
                        basis_file=p.get('basis_file'))

                else:  # Load everything else
                    self.system.add_elementary_particle(
                        symbol, p.get('origin'),
                        basis_name=p.get('basis_name'),
                        basis_file=p.get('basis_file'))

        # set charges
        self.system.set_charges(
            self.data.charges, self.data.scf.get('method', 'rhf') == 'uhf')

    def solve(self):
        """
        Executes the tasks to be performed
        """
        methods = {
           'uhf': napmo.HF,
           'hf':  napmo.HF,
           'rhf': napmo.HF,
           'dft': napmo.DFT
        }

        method = methods.get(self.data.scf.get('method', 'rhf'), None)

        if method is None:
            raise NotImplementedError(self.data.scf.get('method', 'rhf')+" method Not Implemented")
        else:
            self.solver = method(self.system, options=self.data.scf)
            self._energy = self.solver.compute()

    def exec_code(self):
        """
        Executes the source code from the block ``code`` in the input file.
        """
        print('\n--------------------------------------------------')

        g = {'plt': plt, 'print:': print}

        l = {"Energy": self._energy,
             "D": [psi.D for psi in self.solver.PSI],
             "C": [psi.C for psi in self.solver.PSI],
             "E": [psi.O for psi in self.solver.PSI],
             "ND": [psi.D for psi in self.solver.NPSI],
             "NC": [psi.C for psi in self.solver.NPSI],
             "NE": [psi.O for psi in self.solver.NPSI]}

        if self.data.code != '':
            print('Code results:')
            exec(self.data.code, g, l)

    def show(self):
        print(self.system)
