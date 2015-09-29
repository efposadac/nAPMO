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


def second_der_z(r, rm):
    return (rm**2 * (rm + 3.0 * r))/(2.0 * np.pi * ((rm * r)/(rm + r)**2)**(1.5) * (rm + r)**5)


def first_der_z(r, rm):
    return -(np.sqrt(((rm * r)/(rm + r)**2))/(np.pi * r))


def finite_diff_3_matrix(grid, r, rm, l, h):

    A = np.zeros((grid.n_radial+2, grid.n_radial+2))
    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    for i in range(grid.n_radial):
        dzdr = first_der_z(r[i], rm)
        dzdr *= dzdr
        d2zdr2 = second_der_z(r[i], rm)
        aux0 = l * (l + 1.0) / (r[i] * r[i])
        aux1 = dzdr / (h * h)
        aux2 = d2zdr2 / (2.0 * h)

        A[i+1, i] = aux1 - aux2
        A[i+1, i+1] = -(aux0 + (2.0 * aux1))
        A[i+1, i+2] = aux1 + aux2

    return A


def finite_diff_3_rho(grid, r, p_lm, u_00, lindex):
    # Boundary conditions
    p = np.zeros(grid.n_radial+2)
    p[0] = u_00
    for i in range(grid.n_radial):
        p[i+1] = -r[i] * p_lm[i, lindex] * np.pi * 4.0

    return p


def build_potential(grid, U_lm, r, lmax):
    # build the potential
    idx = 0
    V_ab = np.zeros(grid.size)
    for i in range(grid.n_radial):
        aux = 1.0/r[i]
        for j in range(grid.n_angular):
            lindex = 0
            for l in range(lmax+1):
                for m in range(-l, l+1):
                    V_ab[idx] += U_lm[i, lindex] * Y(l, m, grid.angular_theta[j], grid.angular_phi[j]) * aux
                    lindex += 1
            idx += 1

    return V_ab


def coulomb_potential(a, b, grid, molecule, lmax):

    # calculate p_lm
    def p_ab(coord):
        return a.compute(coord) * b.compute(coord)

    rm = molecule.get('atoms')[-1].get('atomic_radii_2')
    r = np.array([grid.radial_abscissas[i] for i in range(grid.n_radial)]) * rm
    p_lm = rho_lm(p_ab, r, angularPoints, lmax)

    # Calculate boundaries
    pi_4 = 4.0 * np.pi
    q_n = grid.integrate(molecule, p_ab)  # integral single center source density
    u_00 = np.sqrt(pi_4 * q_n)
    # print("q_n", q_n, u_00)

    # Solve U_lm
    lsize = (lmax + 1)**2
    U_lm = np.zeros((grid.n_radial, lsize))

    # Obtain z points.
    nradial_2 = grid.n_radial + 2
    z = chebgauss_z(nradial_2)
    h = z[0]

    lindex = 0
    for l in range(lmax+1):
        # Build A matrix
        A = finite_diff_3_matrix(grid, r, rm, l, h)
        # print_matrix(A, nradial_2)
        for m in range(-l, l+1):
            # Build rho
            p = finite_diff_3_rho(grid, r, p_lm, u_00, lindex)
            # Solve the linear problem (get x) (Ax = p)
            x = np.linalg.solve(A, p)
            # Store results
            U_lm[:, lindex] = x[1:-1]

            lindex += 1

    # make interpolation (pending)

    # build the potential
    V_ab = build_potential(grid, U_lm, r, lmax)

    # Some plotting
    # recovered = recover_rho(grid, p_lm, lmax)
    # plt.plot(r, p_ab(grid.xyz)[::angularPoints], 'r-', label='Patron')
    # plt.plot(r, recovered[::angularPoints], 'gx-', label='Recovered')
    # plt.plot(r, V_ab[::grid.n_angular], 'b-', label='potential')
    # plt.legend()
    # plt.xlim([0, 10])
    # plt.show()

    return V_ab


def coulomb_integral(c, d, V_ab, grid, molecule):
    particle = molecule.get('atoms')[-1]
    rm = particle.get('atomic_radii_2')

    r = np.zeros([3], dtype=np.float64)
    integral = 0.0
    for j in range(grid.size):
        r = grid.xyz[j, :]
        p = grid.weight(r, 0, molecule.get('atoms'))
        aux = r - particle.get("origin")
        F = c.compute(r) * d.compute(r) * V_ab[j]
        integral += aux.dot(aux) * p * grid.w[j] * rm * F

    return integral * 4.0 * np.pi


def coulomb_integrals(molecule, radialPoints, angularPoints):
    # For one atom
    lmax = int(lebedev_get_order(angularPoints)/2)

    # Initializing data

    grid = BeckeGrid(radialPoints, angularPoints)
    # print(grid.integrate(molecule, rho_r, args=(molecule,)))

    # Scale and move grids (output in cartesian)
    grid.move(scaling_factor=molecule.get('atoms')[-1].get('atomic_radii_2'))

    # Start integrals calculation
    basis = molecule.get_basis_set('e-')
    length = basis.get('length')
    values = []
    for a in range(length):
        n = a
        for b in range(a, length):
            u = b
            # Calculate V_{ab}
            V_ab = coulomb_potential(basis.get('function')[a], basis.get('function')[b], grid, molecule, lmax)
            for c in range(n, length):
                for d in range(u, length):
                    # Calculate int V_{ab} \phi_c \phi_d dr^3
                    val = coulomb_integral(basis.get('function')[c], basis.get('function')[d], V_ab, grid, molecule)
                    values.append(val)
                    print(a, b, c, d, val)
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
    radialList = [1000]  # [i for i in range(10, 5000, 100)]

    integrals = []
    for radialPoints in radialList:
        integrals.append(coulomb_integrals(molecule, radialPoints, angularPoints))

    # plt.plot(radialList, [0.79788456080286518]*len(radialList))
    # plt.plot(radialList, integrals)
    # plt.show()
