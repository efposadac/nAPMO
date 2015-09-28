#!/usr/bin/env python3
# file: Coulomb.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division

import numpy as np
import scipy as sci
import scipy.special as sp
import matplotlib.pyplot as plt

from napmo.utilities.angular_quadratures import *
from napmo.utilities.radial_quadratures import *
from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *


def print_matrix(A, n):
    for i in range(n):
        for j in range(n):
            print("%12.5f" % (A[i, j]), end="")
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


def rho_n(n, a, b, grid, molecule):
    P_n = np.empty(grid.size)
    r = np.zeros([3], dtype=np.float64)

    particle = molecule.get('atoms')[n]
    rm = particle.get('atomic_radii_2')

    if grid.spherical:
        grid.convert()

    for i in range(grid.size):
        r[0] = grid.x[i]
        r[1] = grid.y[i]
        r[2] = grid.z[i]

        p = grid.weight(r, n, molecule.get('atoms'))
        aux = r - particle.get("origin")
        P_n[i] = aux.dot(aux) * p * grid.w[i] * rm * a.compute(r) * b.compute(r)

    return P_n


def rho_r(coord, molecule):
        basis = molecule.get_basis_set('e-')
        occupation = molecule.n_occupation('e-')
        bvalue = basis.compute(coord)
        output = 0.0

        for k in range(occupation):
            output += bvalue[k] * bvalue[k]

        return output * 2


def second_der_z(r, rm):
    return (rm**2 * (rm + 3.0 * r))/(2.0 * np.pi * ((rm * r)/(rm + r)**2)**(1.5) * (rm + r)**5)


def first_der_z(r, rm):
    return -(np.sqrt(((rm * r)/(rm + r)**2))/(np.pi * r))


def pot_lm(grids):
    pass


def coulomb_potential(a, b, grids, Y_lm, molecule, lmax):
    # calculate p_lm for first atom
    P_n = rho_n(0, a, b, grids[0], molecule)  # P_n[grid.size]
    P_lm = rho_lm(0, a, b, grids[0], Y_lm[0], lmax, P_n, molecule)  # P_lm[grid.n_radial, lm.size]

    # Calculate boundaries
    q_n = P_n.sum() * 4.0 * np.pi  # integral single center source
    u_00 = np.sqrt(4.0 * np.pi * q_n)

    # print("q_n", q_n, u_00)

    # Solve U_lm
    U_lm = np.zeros([grids[0].n_radial, (lmax + 1)**2])

    rm = molecule.get('atoms')[0].get('atomic_radii_2')

    # Obtain z points.
    z = chebgauss_z(grids[0].n_radial + 2)
    rr = np.cos(np.pi * z)
    rr = rm * (1.0 + rr) / (1.0 - rr)
    h = z[0]

    A = np.zeros([grids[0].n_radial + 2, grids[0].n_radial + 2])
    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    rho = np.zeros(grids[0].n_radial + 2)
    rho[0] = u_00

    cons = 4.0 * np.pi

    lm_index = 0
    for l in range(lmax+1):
        for m in range(-l, l+1):
            for i in range(grids[0].n_radial):
                r = grids[0].radial_abscissas[i] * rm
                s = first_der_z(r, rm)
                s *= s
                t = second_der_z(r, rm)
                aux = l * (l + 1.0) / (r * r)
                aux1 = s / (h * h)
                aux2 = t / (2.0 * h)
                A[i+1, i] = aux1 - aux2
                A[i+1, i+1] = -(aux + (2.0 * aux1))
                A[i+1, i+2] = aux1 + aux2

                rho[i + 1] = r * P_lm[i][lm_index] * cons

            # print_matrix(A, grids[0].n_radial+2)

            x = np.linalg.solve(A, rho)

            plt.plot(z[1:-1], rho[1:-1], label=str(lm_index))

            U_lm[:, lm_index] = x[1:-1]

            lm_index += 1

    # # make interpolation
    plt.legend()
    plt.show()

    # # build the potential
    V_total = np.zeros([grids[0].size])

    grid_idx = 0
    for i in range(grids[0].n_radial):
        r = 1.0 / grids[0].radial_abscissas[i] * rm
        for j in range(grids[0].n_angular):
            for lm_index in range((lmax + 1)**2):
                V_total[grid_idx] += Y_lm[0][grid_idx][lm_index] * U_lm[i][lm_index] * r
            grid_idx += 1

    return V_total


def coulomb_integral(n, V_ab, c, d, grid, molecule):

    r = np.zeros([3], dtype=np.float64)

    particle = molecule.get('atoms')[n]
    rm = particle.get('atomic_radii_2')

    if grid.spherical:
        grid.convert()

    integral = 0.0
    for i in range(grid.size):
        r[0] = grid.x[i]
        r[1] = grid.y[i]
        r[2] = grid.z[i]

        p = grid.weight(r, n, molecule.get('atoms'))
        aux = r - particle.get("origin")
        integral += aux.dot(aux) * p * grid.w[i] * rm * c.compute(r) * d.compute(r) * V_ab[i]

    return integral * 4.0 * np.pi


def coulomb_integrals(molecule, radialPoints, angularPoints):
    # For one atom
    lmax = int(lebedev_get_order(angularPoints)/2)

    # Initializing data

    grid = BeckeGrid(radialPoints, angularPoints)
    print(grid.integrate(molecule, rho_r))

    #     # Scale and move grids (output in cartesian)
    #     particle = molecule.get('atoms')[i]
    #     grids[-1].move(particle.get("origin"), particle.get('atomic_radii_2'))

    #     # Calculate spherical harmonics
    #     grids[-1].convert()
    #     Y_lm.append(Y_lmax(lmax, grids[-1]))  # Y_lmax[grid.size, lm.size]
    #     grids[-1].convert()  # to cartesian again

    # # Start integrals calculation
    # basis = molecule.get_basis_set('e-')
    # length = basis.get('length')

    # for a in range(length):
    #     n = a
    #     for b in range(a, length):
    #         u = b
    #         # Calculate V_{ab}
    #         V_ab = coulomb_potential(basis.get('function')[a], basis.get('function')[b], grids, Y_lm, molecule, lmax)
    #         for c in range(n, length):
    #             for d in range(u, length):
    #                 print(a, b, c, d, coulomb_integral(0, V_ab, basis.get('function')[c], basis.get('function')[d], grids[0], molecule))
    #                 # Calculate int V_{ab} \phi_c \phi_d dr^3
    #             u = c + 1


if __name__ == '__main__':
    # Grid definition
    angularPoints = 110
    radialPoints = 30

    # Molecule definition
    basis_file = os.path.join(os.path.dirname(__file__), "TEST.json")
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.3704240745], basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # print(molecule.n_occupation('e-'))
    coulomb_integrals(molecule, radialPoints, angularPoints)
