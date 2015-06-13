from __future__ import division

import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

import os
import sys
import yep

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from interfaces.molecular_system import *
from interfaces.becke_grid import *

particle = "N"
basis_name = "STO-3G"
basis_file = "basis.json"
basis_kind = "GTO"

# Molecule definition
molecule = MolecularSystem()

molecule.add_atom(
                particle, [0.000000, 0.000000, 0.70997005],
                basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
                )

molecule.add_atom(
                particle, [0.000000, 0.000000, -0.70997005],
                basis_kind=basis_kind, basis_name=basis_name, basis_file=basis_file
                )

molecule.show()

# Grid definition
angularPoints = 194
radialPoints = 100

grid = BeckeGrid(radialPoints, angularPoints)
atoms = molecule.get('atoms')

basis = deepcopy(atoms[0].get('basis'))
for i in range(1, len(atoms)):
    basis += atoms[i].get('basis')

# Calculating exact value
exact = molecule.n_particles('e-')

# Build Density Matrix
P = np.array(np.loadtxt('/home/fernando/Downloads/V.35/TESTS/RUN2/dens.dat'), order='F', dtype=np.float64)
P = P.reshape(basis.get('length'), basis.get('length'), order='C')

# Functional definition


def rho(coord, particle_stack, P=P, basis=basis):
    bvalue = basis.compute(coord)
    output = bvalue.dot(P.dot(bvalue))
    return output

# Calculate integral
integral = grid.integrate(atoms, rho)


print(radialPoints, angularPoints, "\int rho P: ", integral, "Error: ", exact - integral)
