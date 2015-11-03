# file: test_multicenter_integrator.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import numpy as np
from copy import deepcopy
import os

from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *


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

    results = np.array([1.9999999595998637, 6.000002321070787, 7.999997958358362,
                        9.999998503574844, 11.999997997974441, 13.999997970177285, 16.00000891580461])
    count = 0

    for (element, distance) in zip(elements, distances):

        print(element)

        # Molecule definition
        molecule = MolecularSystem()
        molecule.add_atom(element, [0.0, 0.0, distance / 2.0],
                          basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)
        molecule.add_atom(element, [0.0, 0.0, -distance / 2.0],
                          basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)

        # Get the stack of atoms.
        atoms = molecule.get('atoms')

        # Build the C interface.
        system = CBinding(atoms)

        # Get the density matrix (from a previous calculation)
        P = np.array(np.loadtxt(element + '_dens.dat'),
                     order='F', dtype=np.float64)
        os.system('cp ' + element + '_dens.dat data.dens')

        # Functional definition (for Python)
        def rho(coord, molecule, P=P):
            basis = molecule.get_basis_set('e-')
            occupation = molecule.n_occupation('e-')
            bvalue = np.array(basis.compute(coord))
            output = np.empty(coord.shape[0])
            for i in range(coord.shape[0]):
                output[i] = bvalue[:, i].dot(P.dot(bvalue[:, i]))
            return output

        # Calculate integral (Python Code)
        integral_p = grid.integrate(molecule, rho, args=(molecule,))

        # Calculate integral (C Code)
        integral_c = grid.integrate_c(system)

        np.testing.assert_allclose(integral_p, results[count])
        np.testing.assert_allclose(integral_c, results[count])

        print(integral_p, integral_c)

        # Delete temporary files.
        os.system('rm data.dens')

        count += 1

    grid.free()

# test_multicenter_integrator()
