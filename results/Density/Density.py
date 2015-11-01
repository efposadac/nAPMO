#!/usr/bin/env python3
# file: Density.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division

import numpy as np
from copy import deepcopy

import os
import time

from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *

"""
Calcualtion of :math:`\int \\rho(\\bf r)` for several diatomic molecules.
"""


def init_system(element, distance, basis_kind, basis_file):
    # Molecule definition
    molecule = MolecularSystem()
    molecule.add_atom(element, [0.000000, 0.000000, distance/2.0], basis_kind=basis_kind, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance/2.0], basis_kind=basis_kind, basis_file=basis_file)
    # molecule.show()

    # Get the stack of atoms.
    atoms = molecule.get('atoms')

    # Build the C interface.
    system = CBinding(atoms)

    # Calculating exact value (Total number of electrons in the system)
    exact = molecule.n_particles('e-')

    # Get the density matrix (from a previous calculation)
    file_dens = os.path.join(os.path.dirname(__file__), element+'_dens.dat')
    P = np.array(np.loadtxt(file_dens), order='F', dtype=np.float64)
    os.system('cp '+file_dens+' data.dens')

    # Functional definition (for Python)
    def rho(coord, molecule, P=P):
        basis = molecule.get_basis_set('e-')
        occupation = molecule.n_occupation('e-')
        bvalue = np.array(basis.compute(coord)).flatten()
        output = bvalue.dot(P.dot(bvalue))
        return output

    return molecule, system, exact, rho


if __name__ == '__main__':

    # Grid definition
    angularPoints = 194
    radialPoints = 100
    grid = BeckeGrid(radialPoints, angularPoints)
    grid.show()

    # Test for diatomic molecules for the following elements:
    elements = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O']
    distances = [0.742, 2.623, 2.427, 1.586, 1.268, 1.098, 1.206]
    basis_name = "STO-3G"
    basis_file = os.path.join(os.path.dirname(__file__), "STO-3G.json")
    basis_kind = "GTO"

    # Header for results.
    print("System Int C         Int Py        Error          Time Py       Time C")

    for (element, distance) in zip(elements, distances):
        molecule, system, exact, rho = init_system(element, distance, basis_kind, basis_file)

        # Calculate integral (Python Code)
        start_time = time.time()
        integral_p = grid.integrate(molecule, rho, args=(molecule,))
        elapsed_time_p = time.time() - start_time

        # Calculate integral (C Code)
        start_time = time.time()
        integral_c = grid.integrate_c(system)
        elapsed_time_c = time.time() - start_time

        # Print the results.
        print("%4s %12.8f  %12.8f  %12.8f  %12.7f  %12.7f" % (
            element+str(2), integral_c, integral_p, np.abs(exact - integral_p), elapsed_time_p, elapsed_time_c))

        # Delete temporary files.
        os.system('rm data.dens')

    grid.free()
