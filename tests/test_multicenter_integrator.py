# file: test_multicenter_integrator.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
from copy import deepcopy
import os
import sys

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.molecular_system import *
from interfaces.becke_grid import *


def test_multicenter_integrator():

    # Grid definition
    angularPoints = 194
    radialPoints = 100
    grid = BeckeGrid(radialPoints, angularPoints)
    grid.show()

    # Test for diatomic molecules of:
    elements = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O']
    distances = [0.742, 2.623, 2.427, 1.586, 1.268, 1.098, 1.206]
    basis_name = "STO-3G"
    basis_file = "STO-3G.json"
    basis_kind = "GTO"

    results = np.array([2.00010688, 6.00025286, 8.00023682, 10.00049433, 12.00127662, 14.00085234, 16.00068683])
    count = 0

    for (element, distance) in zip(elements, distances):

        print(element)

        # Molecule definition
        molecule = MolecularSystem()
        molecule.add_atom(element, [0.0, 0.0, distance/2.0], basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)
        molecule.add_atom(element, [0.0, 0.0, -distance/2.0], basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)

        # Get the stack of atoms.
        atoms = molecule.get('atoms')

        # Build the C interface.
        system = CBinding(atoms)

        # Combine basis-set objects in one.
        basis = deepcopy(atoms[0].get('basis'))
        for i in range(1, len(atoms)):
            basis += atoms[i].get('basis')

        # Get the density matrix (from a previous calculation)
        P = np.array(np.loadtxt(element+'_dens.dat'), order='F', dtype=np.float64)
        os.system('cp '+element+'_dens.dat data.dens')

        # Functional definition (for Python)

        def rho(coord, particle_stack, P=P, basis=basis):
            bvalue = basis.compute(coord)
            output = bvalue.dot(P.dot(bvalue))
            return output

        # Calculate integral (Python Code)
        integral_p = grid.integrate(atoms, rho)

        # Calculate integral (C Code)
        integral_c = grid.integrate_c(system)

        np.testing.assert_allclose(integral_p, results[count])
        np.testing.assert_allclose(integral_c, results[count])

        # Delete temporary files.
        os.system('rm data.dens')

        count += 1

    grid.free()

# test_multicenter_integrator()
