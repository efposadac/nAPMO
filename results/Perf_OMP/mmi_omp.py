#!/usr/bin/env python3
# file: mmi_omp.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke_grid import BeckeGrid


def init_system(element, distance, basis_kind, basis_file):
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom(element, [0.000000, 0.000000, distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)
    # molecule.show()

    # Get the stack of atoms.
    atoms = molecule.get('atoms')

    # Build the C interface.
    system = CBinding(atoms)

    # Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')

    return atoms, system, exact


if __name__ == '__main__':

    # Grid definition
    angularPoints = 5810
    radialPoints = 1000
    grid = BeckeGrid(radialPoints, angularPoints)

    atoms, system, exact = init_system('C', 1.268, 'GTO', 'STO-3G.json')

    # Calculate integral (C Code)
    integral = grid.integrate_c(system)

    grid.free()
