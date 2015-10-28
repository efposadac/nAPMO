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
import scipy as sci
import scipy.special as sp
import matplotlib.pyplot as plt
import time

from napmo.utilities.angular_quadratures import *
from napmo.utilities.radial_quadratures import *
from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *

from horton import *


def print_matrix(A, n):
    for i in range(n):
        for j in range(n):
            print("%12.5f" % (A[i, j]), end="")
        print("")
    print("")


def real_spherical_harmonics(m, l, theta, phi):
    if m == 0:
        aux = sp.sph_harm(m, l, theta, phi)
    elif m < 0:
        aux_a = sp.sph_harm(-m, l, theta, phi)
        aux_b = sp.sph_harm(m, l, theta, phi)
        aux = 1.0j * 0.70710678118654757 * (aux_b + aux_a)
    else:
        aux_a = sp.sph_harm(m, l, theta, phi)
        aux_b = sp.sph_harm(-m, l, theta, phi)
        aux = 0.70710678118654757 * (aux_b - aux_a)
    return aux.real


def Y(l, m, theta, phi):
    return real_spherical_harmonics(m, l, phi, theta)


# Define integrand P(r, theta, phi)
def int_l(theta, phi, func, r, l, m):

    # Convert to cartesian
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    coord = np.dstack((x, y, z)).reshape(len(x), 3)

    return func(coord) * Y(l, m, theta, phi)


def rho_lm(func, rad, n, lmax):
    # Size of expansion
    lsize = (lmax + 1) * (lmax + 1)

    # Start calculation
    p_lm = np.zeros((len(rad), lsize))

    for r in range(len(rad)):
        lindex = 0
        for l in range(lmax+1):
            for m in range(-l, l+1):
                # integrate
                p_lm[r, lindex] = lebedev_integrate(int_l, n, args=(func, rad[r], l, m))
                lindex += 1

    return p_lm


def recover_rho(grid, p_lm, lmax):
    idx = 0
    rho = np.zeros(grid.size)
    for r in range(grid.n_radial):
        for a in range(grid.n_angular):
            lindex = 0
            for l in range(lmax+1):
                for m in range(-l, l+1):
                    rho[idx] += p_lm[r, lindex] * Y(l, m, grid.angular_theta[a], grid.angular_phi[a])
                    lindex += 1
            idx += 1

    return rho


def rho_r(coord, molecule):
        basis = molecule.get_basis_set('e-')
        occupation = molecule.n_occupation('e-')
        bvalue = basis.compute(coord)
        output = 0.0

        for k in range(occupation):
            output += bvalue[k] * bvalue[k]

        return output * 2


def r_to_z(r, rm):
    return np.arccos((r - rm)/(r + rm))/np.pi


def second_der_z(r, rm):
    return (rm**2 * (rm + (3.0 * r)))/(2.0 * np.pi * (((rm * r)/(rm + r)**2)**(1.5)) * (rm + r)**5)


def first_der_z(r, rm):
    return -(np.sqrt(((rm * r)/(rm + r)**2))/(np.pi * r))


def finite_diff_7m(grid, r, rm, l, h):
    aux_h2 = 1.0 / (h * h)
    aux_h = 1.0 / h

    A = np.zeros((grid.n_radial+2, grid.n_radial+2))

    for i in range(3):
        A[i, i] = 1.0
        A[-1-i, -1-i] = 1.0

    # i = 1
    f_der_coeff = np.array(
        [-0.166666666666666, -1.283333333333333, 2.5, -1.666666666666666, 0.833333333333333, -0.2499999999999999, 0.03333333333333333, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.7, -0.3888888888888888, -2.7, 4.75, -3.722222222222222, 1.8, -0.5, 0.06111111111111111],
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
        [-0.033333333333333215, 0.24999999999999956, -0.8333333333333326, 1.6666666666666659, -2.499999999999999, 1.283333333333333, 0.16666666666666666, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.06111111111111, -0.49999999999999, 1.7999999999999, -3.7222222222222, 4.74999999999999, -2.69999999999999, -0.38888888888888, 0.7],
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
        A[-2, -i-1] = s_der_coeff[-i-1] + f_der_coeff[-i-1]

    A[-2, -2] += -aux0

    # forward i = 2

    f_der_coeff = np.array(
        [0.033333333333333, -0.4, -0.583333333333333, 1.3333333333333333, -0.5, 0.13333333333333333, -0.016666666666666666, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [-0.06111111111111111, 1.188888888888888, -2.1, 0.7222222222222222, 0.4722222222222222, -0.3, 0.0888888888888888, -0.01111111111111111],
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
        [0.0166666666666666, -0.133333333333333, 0.499999999999999, -1.333333333333333, 0.5833333333333333, 0.3999999999999999, -0.03333333333333333, 0.0],
        dtype=np.float64)

    s_der_coeff = np.array(
        [-0.011111111111111, 0.088888888888888, -0.299999999999999, 0.47222222222222, 0.722222222222222, -2.099999999999999, 1.188888888888888, -0.06111111111111111],
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
        A[-3, -i-1] = s_der_coeff[-i-1] + f_der_coeff[-i-1]

    A[-3, -3] += -aux0

    # 7 points

    f_der_coeff = np.array(
        [-0.016666666666666607, 0.1499999999999999, -0.75, 1.1102230246251565e-16, 0.7499999999999999, -0.15, 0.016666666666666666],
        dtype=np.float64)

    s_der_coeff = np.array(
        [0.011111111111111072, -0.1499999999999997, 1.4999999999999993, -2.7222222222222214, 1.4999999999999998, -0.14999999999999997, 0.011111111111111112],
        dtype=np.float64)

    for i in range(grid.n_radial-4):
        dzdr = first_der_z(r[i+2], rm)
        dzdr *= dzdr
        d2zdr2 = second_der_z(r[i+2], rm)

        aux0 = l * (l + 1.0) / (r[i+2] * r[i+2])
        aux1 = dzdr * aux_h2
        aux2 = d2zdr2 * aux_h

        s_der_coeff_aux = s_der_coeff * aux1
        f_der_coeff_aux = f_der_coeff * aux2

        for j in range(len(f_der_coeff_aux)):
            A[i+3, i+j] = s_der_coeff_aux[j] + f_der_coeff_aux[j]

        A[i+3, i+3] += -aux0

    return A


def finite_diff_7r(grid, r, p_lm, u_00, l, lindex):
    # Boundary conditions
    p = np.zeros(grid.n_radial+2)
    for i in range(grid.n_radial):
        p[i+1] = -r[i] * p_lm[i, lindex] * np.pi * 4.0

    if l == 0:
        p[0] = u_00

    return p


def finite_diff_3m(grid, r, rm, l, h):
    aux_h2 = 1.0 / (h * h)
    aux_h = 1.0 / h

    A = np.zeros((grid.n_radial+2, grid.n_radial+2))
    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    f_der_coeff = np.array(
        [-0.5, 0.0, 0.5],
        dtype=np.float64)

    s_der_coeff = np.array(
        [1.0, -2.0, 1.0],
        dtype=np.float64)

    for i in range(grid.n_radial):
        dzdr = first_der_z(r[i], rm)
        dzdr *= dzdr
        d2zdr2 = second_der_z(r[i], rm)

        aux0 = l * (l + 1.0) / (r[i] * r[i])
        aux1 = dzdr * aux_h2
        aux2 = d2zdr2 * aux_h

        s_der_coeff_aux = s_der_coeff * aux1
        f_der_coeff_aux = f_der_coeff * aux2

        for j in range(len(f_der_coeff_aux)):
            A[i+1, i+j] = s_der_coeff_aux[j] + f_der_coeff_aux[j]

        A[i+1, i+1] += -aux0

    return A


def finite_diff_3r(grid, r, p_lm, u_00, l, lindex):
    # Boundary conditions
    p = np.zeros(grid.n_radial+2)
    for i in range(grid.n_radial):
        p[i+1] = -r[i] * p_lm[i, lindex] * np.pi * 4.0

    p[0] = u_00

    return p


def build_potential(grid, U_lm, r, lmax):
    # build the potential
    idx = 0
    V_ab = np.zeros(grid.size)
    for i in range(grid.n_radial):
        aux = 1 / r[i]
        for j in range(grid.n_angular):
            lindex = 0
            for l in range(lmax+1):
                for m in range(-l, l+1):
                    V_ab[idx] += U_lm[i, lindex] * Y(l, m, grid.angular_theta[j], grid.angular_phi[j]) * aux
                    lindex += 1
            idx += 1

    return V_ab


def solve_poisson(p_lm, lmax, molecule, grid, dens, fd=3):
    rm = molecule.get('atoms')[-1].get('atomic_radii_2')
    r = np.array([grid.radial_abscissas[i] for i in range(grid.n_radial)]) * rm

    # Calculate boundaries
    q_n = grid.integrate(molecule, dens)  # integral single center source density

    u_00 = 0.0
    if np.abs(q_n) > 1.0e-16:
        u_00 = np.sqrt(4.0 * np.pi * q_n)

    # Solve U_lm
    lsize = (lmax + 1)**2
    U_lm = np.zeros((grid.n_radial, lsize))

    # Obtain z points.
    z = r_to_z(r, rm)
    h = z[0]

    # Solve U_lm
    lsize = (lmax + 1)**2
    U_lm = np.zeros((grid.n_radial, lsize))

    # Obtain z points.
    z = r_to_z(r, rm)
    h = z[0]

    lindex = 0
    for l in range(lmax+1):
        # Build A matrix
        if fd == 3:
            A = finite_diff_3m(grid, r, rm, l, h)
        elif fd == 7:
            A = finite_diff_7m(grid, r, rm, l, h)

        # print_matrix(A, grid.n_radial+2)
        for m in range(-l, l+1):
            # Build rho
            if fd == 3:
                p = finite_diff_3r(grid, r, p_lm, u_00, l, lindex)
            elif fd == 7:
                p = finite_diff_7r(grid, r, p_lm, u_00, l, lindex)

            # Solve the linear problem (get x) (Ax = p)
            x = np.linalg.solve(A, p)

            # Store results
            U_lm[:, lindex] = x[1:-1]

            lindex += 1

    return U_lm


def coulomb_potential(a, b, grid, molecule, lmax, fd):

    def p_ab(coord):
        return a.compute(coord) * b.compute(coord)

    # calculate p_lm
    rm = molecule.get('atoms')[-1].get('atomic_radii_2')
    rad = np.array([grid.radial_abscissas[i] for i in range(grid.n_radial)]) * rm
    p_lm = rho_lm(p_ab, rad, angularPoints, lmax)

    # Solve poisson for p_lm
    U_lm = solve_poisson(p_lm, lmax, molecule, grid, p_ab, fd=fd)

    # Build potential
    V_ab = build_potential(grid, U_lm, rad, lmax)

    return V_ab


def coulomb_integral(c, d, V_ab, grid, molecule):
    particle = molecule.get('atoms')[-1]
    rm = particle.get('atomic_radii_2')

    if grid.spherical:
        return None

    r = np.zeros([3], dtype=np.float64)
    integral = 0.0
    for j in range(grid.size):
        r = grid.xyz[j, :]
        p = grid.weight(r, 0, molecule.get('atoms'))
        aux = r - particle.get("origin")
        F = c.compute(r) * d.compute(r) * V_ab[j]
        integral += aux.dot(aux) * p * grid.w[j] * rm * F

    return integral * 4.0 * np.pi


def coulomb_integrals(molecule, radialPoints, angularPoints, fd):
    # For one atom
    lmax = int(lebedev_get_order(angularPoints)/2)

    # Initializing data
    grid = BeckeGrid(radialPoints, angularPoints)

    # Scale and move grids (output in cartesian)
    grid.move(scaling_factor=molecule.get('atoms')[-1].get('atomic_radii_2'))

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
            V_ab = coulomb_potential(basis.get('function')[a], basis.get('function')[b], grid, molecule, lmax, fd)
            for c in range(n, length):
                for d in range(u, length):
                    # Calculate int V_{ab} \phi_c \phi_d dr^3
                    val = coulomb_integral(basis.get('function')[c], basis.get('function')[d], V_ab, grid, molecule)
                    print(a+1, b+1, c+1, d+1, val)
                u = c + 1

    return val

    grid.free()


if __name__ == '__main__':

    # Molecule definition
    basis_file = os.path.join(os.path.dirname(__file__), "TEST.json")
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0], basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # Grid definition
    angularPoints = 14
    # radialList = [10]
    # radialList = [i for i in range(6, 1000, 10)]
    radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    # fd_alg = [765]
    fd_alg = [3, 7]
    for i in fd_alg:
        integrals = []
        times = []
        for radialPoints in radialList:
            print(i, radialPoints)
            start_time = time.time()
            integrals.append(coulomb_integrals(molecule, radialPoints, angularPoints, i))
            times.append(time.time() - start_time)

        plt.yscale('log')
        plt.xlabel('Radial Points')
        plt.ylabel('Error (log)')
        plt.plot(radialList, np.abs(np.array(integrals)-0.79788456080286518), label=str(i))

    plt.legend(loc=2)
    plt.show()
