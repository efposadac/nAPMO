#!/usr/bin/env python2

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

from Coulomb import *

if __name__ == '__main__':

    basis_file = "/home/fernando/PhD/Develop/nAPMO/results/Coulomb/TEST.json"
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0], basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    rm = molecule.get('atoms')[-1].get('atomic_radii_2')

    # Grid definition
    angularPoints = 14
    radialList = [100]
    # radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    h_ints = []
    n_ints = []
    n2_ints = []
    for radialPoints in radialList:
        lmax = int(lebedev_get_order(angularPoints)/2)
        lsize = (lmax + 1) * (lmax + 1)

        # Initializing data

        grid = BeckeGrid(radialPoints, angularPoints)
        grid.move()
        r = grid.radial_abscissas.copy()
        # r *= rm
        print(rm)

        print(np.allclose(grid.radial_abscissas, r))

        h_grid_spec = AtomicGridSpec('power:'+str(r.min())+':'+str(r.max())+':'+str(radialPoints)+':'+str(angularPoints))
        molgrid = BeckeMolGrid(np.array([[0.0, 0.0, 0.0]]), np.array([2]), np.array([2.0]), agspec=h_grid_spec, random_rotate=False, mode='keep')
        atgrid = molgrid.subgrids[-1]
        rad = atgrid.rgrid.radii

        basis = molecule.get('atoms')[-1].get('basis')
        a = basis.get('function')[0]
        b = basis.get('function')[0]

        a.show()
        b.show()

        def p_ab(coord):
            return a.compute(coord) * b.compute(coord)

        # #############################
        # ### Spherical expansion test
        # #############################

        # npamo
        p_lm = rho_lm(p_ab, rad, angularPoints, lmax)

        # horton
        becke_weights = molgrid.becke_weights

        rho = np.empty(molgrid.size)
        for i in range(molgrid.size):
            rho[i] = p_ab(molgrid.points[i, :])

        density_decomposition = atgrid.get_spherical_decomposition(rho, becke_weights, lmax=lmax)

        p_lm_h = np.empty(p_lm.shape)
        for i in range(lsize):
            p_lm_h[:, i] = density_decomposition[i].y[:]

        # compare horton - napmo
        print('spherical decomposition:', np.allclose(p_lm, p_lm_h))

        # Build rho again:

        # napmo
        rho_napmo = build_potential(grid, p_lm, np.ones(grid.size), lmax)
        print('rebuild napmo:', np.allclose(rho_napmo, rho))

        # horton
        rho_horton = molgrid.zeros()
        molgrid.eval_decomposition(density_decomposition, molgrid.centers[-1], rho_horton)
        print('rebuild horton:', np.allclose(rho_horton, rho))

        # #############################
        # ### Integration test
        # #############################

        # Solve Poisson
        # #############

        # horton
        U_lm_h = solve_poisson_becke(density_decomposition)
        # U_lm_n = np.empty(p_lm.shape)
        # for i in range(lsize):
        #     U_lm_n[:, i] = U_lm_h[i].y[:]

        # napmo
        p_lm = rho_lm(p_ab, r, angularPoints, lmax)
        U_lm_n = solve_poisson(p_lm, lmax, molecule, grid, p_ab, fd=7)

        # for i in range(lsize):
        #     plt.plot(r[::-1], U_lm_n[:, i], '-x', label='napmo'+str(i))
        #     plt.plot(rad, U_lm_h[i].y, '-x', label='horton'+str(i))
        # plt.legend()
        # plt.show()

        # Build the potential and calculate integral
        # ##########################################

        # napmo
        V_ab_n = build_potential(grid, U_lm_n, r, lmax)
        val_n = coulomb_integral(a, b, V_ab_n, grid, molecule)
        n_ints.append(val_n)

        # horton
        V_ab_h = molgrid.zeros()
        molgrid.eval_decomposition(U_lm_h, molgrid.centers[-1], V_ab_h)
        val_h = atgrid.integrate(rho, V_ab_h)
        h_ints.append(val_h)

        # print(V_ab_n)
        # plt.plot(r, V_ab_n[::angularPoints], '-x', label='napmo')
        # plt.plot(rad, V_ab_h[::angularPoints], '-x', label='horton')
        # plt.legend()
        # plt.show()

        U_lm_n = solve_poisson(p_lm, lmax, molecule, grid, p_ab, fd=3)
        V_ab_n = build_potential(grid, U_lm_n, r, lmax)
        val_n2 = coulomb_integral(a, b, V_ab_n, grid, molecule)
        n2_ints.append(val_n2)

        print('Potential:', np.allclose(V_ab_h, V_ab_n))
        print(rm, val_n2, val_n)

    plt.yscale('log')
    plt.xlabel('Radial Points')
    plt.ylabel('Error (log)')
    # plt.plot(radialList, np.abs(np.array(h_ints)-0.79788456080286518), label='horton')
    plt.plot(radialList, np.abs(np.array(n_ints)-0.79788456080286518), label='napmo-7')
    plt.plot(radialList, np.abs(np.array(n2_ints)-0.79788456080286518), label='napmo-3')
    plt.legend()
    plt.savefig('error.png')
