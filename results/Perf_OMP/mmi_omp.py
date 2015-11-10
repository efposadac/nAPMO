#!/usr/bin/env python3
# file: mmi_omp.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.utilities.density import *
import os


def init_system(element, distance, basis_kind, basis_file):
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom(element, [0.000000, 0.000000, distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance / 2.0],
                      basis_kind=basis_kind, basis_file=basis_file)
    # molecule.show()
    exact = molecule.n_particles('e-')

    return molecule, exact


if __name__ == '__main__':

    # Grid definition
    angularPoints = 5810
    radialPoints = 1000

    molecule, exact = init_system('C', 1.268, 'GTO', 'STO-3G.json')
    grid = BeckeGrid(molecule, radialPoints, angularPoints)

    basis = molecule.get_basis_as_cstruct('e-')

    # Calculate integral
    file_dens = os.path.join(os.path.dirname(__file__), 'data.dens')
    f = density_full_from_matrix_gto(file_dens, basis, grid.points)
    integral = grid.integrate(f)
