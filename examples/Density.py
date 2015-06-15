from __future__ import division

import numpy as np
from copy import deepcopy

import os
import sys
import time

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.molecular_system import *
from interfaces.becke_grid import *

# Grid definition
angularPoints = 194
radialPoints = 100

grid = BeckeGrid(radialPoints, angularPoints)

# Test for diatomic molecules for the following elements:
elements = ['H', 'Li', 'Be', 'B', 'C', 'N', 'O']
distances = [0.742, 2.623, 2.427, 1.586, 1.268, 1.098, 1.206]
basis_name = "STO-3G"
basis_file = "STO-3G.json"
basis_kind = "GTO"

for (element, distance) in zip(elements, distances):
    # Molecule definition
    molecule = MolecularSystem()

    molecule.add_atom(element, [0.000000, 0.000000, distance/2.0], basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)
    molecule.add_atom(element, [0.000000, 0.000000, -distance/2.0], basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file)

    # molecule.show()

    atoms = molecule.get('atoms')

    basis = deepcopy(atoms[0].get('basis'))
    for i in range(1, len(atoms)):
        basis += atoms[i].get('basis')

    # Calculating exact value
    exact = molecule.n_particles('e-')

    # Build Density Matrix
    P = np.array(np.loadtxt(element+'_dens.dat'), order='F', dtype=np.float64)

    # Functional definition

    def rho(coord, particle_stack, P=P, basis=basis):
        bvalue = basis.compute(coord)
        output = bvalue.dot(P.dot(bvalue))
        return output

    # Calculate integral
    start_time = time.time()
    integral = grid.integrate(atoms, rho)
    elapsed_time = time.time() - start_time

    np.set_printoptions(precision=3)
    print("NF: %5d %s Int: %8.3f Error: %8.4f Time: %12.7f" % (
        basis.get('t_length'), element+str(2), integral, np.abs(exact - integral), elapsed_time))
