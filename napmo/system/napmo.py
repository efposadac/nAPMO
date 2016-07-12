# file: napmo.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

import napmo


class NAPMO(object):
    """
    nAPMO system manager
    """

    Energy = 0.0

    def __init__(self, data):
        super(NAPMO, self).__init__()

        assert isinstance(data, napmo.InputParser)

        self.data = data

        # Build system
        self.build_system()

        # Show system info
        self.show()

    def build_system(self):
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
        # TODO: Parse method and others (class solver maybe?)
        self.solver = napmo.SCF(self.system, options=self.data.scf)
        return self.solver.compute()

    def show(self):
        print(self.system)
