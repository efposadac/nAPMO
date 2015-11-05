#!/usr/bin/env python
# file: Performance.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

import os
import time

from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke_grid import BeckeGrid


# Molecule definition
def init_system():
    molecule = MolecularSystem()
    molecule.add_atom('H', [0.0, 0.0,  0.371],
                      basis_kind='GTO', basis_file='basis.json')
    molecule.add_atom('H', [0.0, 0.0, -0.371],
                      basis_kind='GTO', basis_file='basis.json')
    molecule.show()

    # Get the stack of atoms.
    atoms = molecule.get('atoms')

    # Build the C interface.
    system = CBinding(atoms)

    # Combine basis-set objects in one.
    basis = deepcopy(atoms[0].get('basis'))
    for i in range(1, len(atoms)):
        basis += atoms[i].get('basis')

    # Get the density matrix (from a previous calculation)
    P = np.array(np.loadtxt('data.dens'), dtype=np.float64)

    # Functional definition (for Python)
    def rho(coord, particle_stack, P=P, basis=basis):
        bvalue = basis.compute(coord)
        output = bvalue.dot(P.dot(bvalue))
        return output

    return atoms, system, rho


# Scaling with respect to the number of radial points.
def radial_scaling():
    atoms, system, rho = init_system()
    # Grid definition
    angularPoints = 194
    radialPoints = [i for i in range(10, 1000, 50)]

    elapsed_time_p = Stack()
    elapsed_time_c = Stack()

    for radialPoint in radialPoints:
        grid = BeckeGrid(radialPoint, angularPoints)

        # Calculate integral (Python Code)
        start_time = time.time()
        integral_p = grid.integrate(atoms, rho)
        elapsed_time_p.push(time.time() - start_time)

        # Calculate integral (C Code)
        start_time = time.time()
        integral_c = grid.integrate_c(system)
        elapsed_time_c.push((time.time() - start_time) * 100)

        # Print the results.
        print("%5d %12.8f  %12.8f  %12.7f  %12.7f" % (
            radialPoint, integral_c, integral_p, elapsed_time_p.peek(), elapsed_time_c.peek()))

        grid.free()

    plt.title('Scaling with the number of radial points')
    plt.xlabel('Number of radial points')
    plt.ylabel('C (cs), Py (s)')
    plt.plot(radialPoints, elapsed_time_c, label='C')
    plt.plot(radialPoints, elapsed_time_p, label='Py')
    plt.legend()
    plt.savefig('radial_points_scaling.png')
    plt.close()


# Scaling with respect to the number of angular points.
def angular_scaling():
    atoms, system, rho = init_system()
    # Grid definition
    angularPoints = [110, 170, 194, 230, 266,
                     302, 350, 434, 590, 770, 974, 1202]
    radialPoints = 100

    elapsed_time_p = Stack()
    elapsed_time_c = Stack()

    for angularPoint in angularPoints:
        grid = BeckeGrid(radialPoints, angularPoint)

        # Calculate integral (Python Code)
        start_time = time.time()
        integral_p = grid.integrate(atoms, rho)
        elapsed_time_p.push(time.time() - start_time)

        # Calculate integral (C Code)
        start_time = time.time()
        integral_c = grid.integrate_c(system)
        elapsed_time_c.push((time.time() - start_time) * 100)

        # Print the results.
        print("%5d %12.8f  %12.8f  %12.7f  %12.7f" % (
            angularPoint, integral_c, integral_p, elapsed_time_p.peek(), elapsed_time_c.peek()))

        grid.free()

    plt.title('Scaling with the number of angular points')
    plt.xlabel('Number of angular points')
    plt.ylabel('C (cs), Py (s)')
    plt.plot(angularPoints, elapsed_time_c, label='C')
    plt.plot(angularPoints, elapsed_time_p, label='Py')
    plt.legend()
    plt.savefig('angular_points_scaling.png')
    plt.close()


if __name__ == '__main__':
    radial_scaling()
    angular_scaling()
