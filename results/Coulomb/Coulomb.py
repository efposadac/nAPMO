#!/usr/bin/env python2
# file: Coulomb.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import time
import os

from napmo.utilities.angular_quadratures import *
from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.grids.poisson_solver import poisson_solver


def finite_diff_7m(grid, r, rm, l, h):
    aux_h2 = 1.0 / (h * h)
    aux_h = 1.0 / h

    A = np.zeros((grid.n_radial + 2, grid.n_radial + 2))

    for i in range(3):
        A[i, i] = 1.0
        A[-1 - i, -1 - i] = 1.0

    # i = 1
    f_der_coeff = np.array(
        [-0.166666666666666, -1.283333333333333, 2.5, -1.666666666666666,
            0.833333333333333, -0.2499999999999999, 0.03333333333333333, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.7, -0.3888888888888888, -2.7, 4.75, -
            3.722222222222222, 1.8, -0.5, 0.06111111111111111],
        dtype=np.float64)

    # i = 1
    dzdr = first_der_z(r[0], rm)
    dzdr *= dzdr
    d2zdr2 = second_der_z(r[0], rm)
    aux0 = l * (l + 1.0) / (r[0] * r[0])
    aux1 = dzdr * aux_h2
    aux2 = d2zdr2 * aux_h

    s_der_coeff *= aux1
    f_der_coeff *= aux2

    for i in range(len(s_der_coeff)):
        A[1, i] = s_der_coeff[i] + f_der_coeff[i]

    A[1, 1] += -aux0

    # i = N-1

    f_der_coeff = np.array(
        [-0.033333333333333215, 0.24999999999999956, -0.8333333333333326,
            1.6666666666666659, -2.499999999999999, 1.283333333333333, 0.16666666666666666, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.06111111111111, -0.49999999999999, 1.7999999999999, -3.7222222222222,
            4.74999999999999, -2.69999999999999, -0.38888888888888, 0.7],
        dtype=np.float64)

    dzdr = first_der_z(r[-1], rm)
    dzdr *= dzdr
    d2zdr2 = second_der_z(r[-1], rm)
    aux0 = l * (l + 1.0) / (r[-1] * r[-1])
    aux1 = dzdr * aux_h2
    aux2 = d2zdr2 * aux_h

    s_der_coeff *= aux1
    f_der_coeff *= aux2

    for i in range(len(s_der_coeff)):
        A[-2, -i - 1] = s_der_coeff[-i - 1] + f_der_coeff[-i - 1]

    A[-2, -2] += -aux0

    # forward i = 2

    f_der_coeff = np.array(
        [0.033333333333333, -0.4, -0.583333333333333, 1.3333333333333333, -
            0.5, 0.13333333333333333, -0.016666666666666666, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [-0.06111111111111111, 1.188888888888888, -2.1, 0.7222222222222222,
            0.4722222222222222, -0.3, 0.0888888888888888, -0.01111111111111111],
        dtype=np.float64)

    dzdr = first_der_z(r[1], rm)
    dzdr *= dzdr
    d2zdr2 = second_der_z(r[1], rm)
    aux0 = l * (l + 1.0) / (r[1] * r[1])
    aux1 = dzdr * aux_h2
    aux2 = d2zdr2 * aux_h

    s_der_coeff *= aux1
    f_der_coeff *= aux2

    for i in range(len(s_der_coeff)):
        A[2, i] = s_der_coeff[i] + f_der_coeff[i]

    A[2, 2] += -aux0

    # i = N - 2

    f_der_coeff = np.array(
        [0.0166666666666666, -0.133333333333333, 0.499999999999999, -1.333333333333333,
            0.5833333333333333, 0.3999999999999999, -0.03333333333333333, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [-0.011111111111111, 0.088888888888888, -0.299999999999999, 0.47222222222222,
            0.722222222222222, -2.099999999999999, 1.188888888888888, -0.06111111111111111],
        dtype=np.float64)

    dzdr = first_der_z(r[-2], rm)
    dzdr *= dzdr
    d2zdr2 = second_der_z(r[-2], rm)
    aux0 = l * (l + 1.0) / (r[-2] * r[-2])
    aux1 = dzdr * aux_h2
    aux2 = d2zdr2 * aux_h

    s_der_coeff *= aux1
    f_der_coeff *= aux2

    for i in range(len(s_der_coeff)):
        A[-3, -i - 1] = s_der_coeff[-i - 1] + f_der_coeff[-i - 1]

    A[-3, -3] += -aux0

    # 7 points

    f_der_coeff = np.array(
        [-0.016666666666666607, 0.1499999999999999, -0.75, 1.1102230246251565e-16,
            0.7499999999999999, -0.15, 0.016666666666666666],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.011111111111111072, -0.1499999999999997, 1.4999999999999993, -2.7222222222222214,
            1.4999999999999998, -0.14999999999999997, 0.011111111111111112],
        dtype=np.float64)

    for i in range(grid.n_radial - 4):
        dzdr = first_der_z(r[i + 2], rm)
        dzdr *= dzdr
        d2zdr2 = second_der_z(r[i + 2], rm)

        aux0 = l * (l + 1.0) / (r[i + 2] * r[i + 2])
        aux1 = dzdr * aux_h2
        aux2 = d2zdr2 * aux_h

        s_der_coeff_aux = s_der_coeff * aux1
        f_der_coeff_aux = f_der_coeff * aux2

        for j in range(len(f_der_coeff_aux)):
            A[i + 3, i + j] = s_der_coeff_aux[j] + f_der_coeff_aux[j]

        A[i + 3, i + 3] += -aux0

    return A


def finite_diff_7r(grid, r, p_lm, u_00, l, lindex):
    # Boundary conditions
    p = np.zeros(grid.n_radial + 2)
    for i in range(grid.n_radial):
        p[i + 1] = -r[i] * p_lm[i, lindex] * np.pi * 4.0

    if l == 0:
        p[0] = u_00

    return p


def coulomb_integrals(molecule, radialPoints, angularPoints, fd):
    # For one atom
    lmax = int(lebedev_get_order(angularPoints) / 2)

    # Initializing data
    grid = BeckeGrid(radialPoints, angularPoints)

    # Scale and move grids (output in cartesian)
    grid.move()

    # Start integrals calculation
    basis = molecule.get_basis_set('e-')
    # basis.show()
    length = basis.get('length')
    values = []
    for a in range(length):
        n = a
        for b in range(a, length):
            u = b
            # Calculate V_{ab}
            V_ab = coulomb_potential(basis.get('function')[a], basis.get(
                'function')[b], grid, molecule, lmax, fd)
            for c in range(n, length):
                for d in range(u, length):
                    # Calculate int V_{ab} \phi_c \phi_d dr^3
                    val = coulomb_integral(basis.get('function')[c], basis.get(
                        'function')[d], V_ab, grid, molecule)
                    print(a + 1, b + 1, c + 1, d + 1, val)
                u = c + 1

    return val

    grid.free()


def rho_ab(a, b, coord):
    return a.compute(coord) * b.compute(coord)

if __name__ == '__main__':

    # Molecule definition
    basis_file = os.path.join(os.path.dirname(__file__), "TEST.json")
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0],
                      basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # basis set
    basis = molecule.get('atoms')[-1].get('basis')
    a = basis.get('function')[0]
    b = basis.get('function')[0]

    # Build grid
    angularPoints = 110
    radialPoints = 55

    grid = BeckeGrid(molecule, radialPoints, angularPoints)

    # Calculate potential
    lmax = int(lebedev_get_order(angularPoints) / 2)
    rho = rho_ab(a, b, grid.points)

    U = poisson_solver(grid, rho, lmax)

    # TODO: build molecular potential
    for i in range(grid.atgrids[-1].radial_grid.size):
        U[i, :] /= grid.atgrids[-1].radial_grid.points[i]

    v = grid.atgrids[-1].evaluate_expansion(lmax, U)

    integral = grid.integrate(rho * v)

    print(integral)
