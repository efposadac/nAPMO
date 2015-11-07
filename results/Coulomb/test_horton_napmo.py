#!/usr/bin/env python2

from __future__ import division
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
import time

from napmo.utilities.angular_quadratures import lebedev_get_order
from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.grids.poisson_solver import poisson_solver

from horton import *

if __name__ == '__main__':

    basis_file = "/home/fernando/PhD/Develop/nAPMO/results/Coulomb/TEST.json"
    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.0],
                      basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    # Grid definition
    angularPoints = 110
    # radialList = [100]
    radialList = [6, 14, 110, 170, 194, 230, 266, 302, 350, 434]

    h_ints = []
    n_ints = []
    n2_ints = []
    for radialPoints in radialList:
        print(radialPoints)
        lmax = int(lebedev_get_order(angularPoints) / 2)
        lsize = (lmax + 1) ** 2

        # Grids
        grid = BeckeGrid(molecule, radialPoints, angularPoints)

        h_grid_spec = AtomicGridSpec('power:' + str(grid.atgrids[-1].radial_grid.points.min()) + ':' + str(
            grid.atgrids[-1].radial_grid.points.max()) + ':' + str(radialPoints) + ':' + str(angularPoints))

        molgrid = BeckeMolGrid(np.array([[0.0, 0.0, 0.0]]), np.array([2]), np.array(
            [2.0]), agspec=h_grid_spec, random_rotate=False, mode='keep')

        atgrid = molgrid.subgrids[-1]

        # Basis and functional
        basis = molecule.get('atoms')[-1].get('basis')
        a = basis.get('function')[0]
        b = basis.get('function')[0]

        def p_ab(coord):
            return a.compute(coord) * b.compute(coord)

        #############################
        # Spherical expansion test
        #############################

        rho = np.empty(molgrid.size)
        rho = p_ab(molgrid.points)

        # napmo
        p_lm_n = grid.atgrids[-1].spherical_expansion(lmax, rho)

        # horton
        becke_weights = molgrid.becke_weights
        density_decomposition = atgrid.get_spherical_decomposition(
            rho, becke_weights, lmax=lmax)

        p_lm_h = np.empty(p_lm_n.shape)
        for i in range(lsize):
            p_lm_h[:, i] = density_decomposition[i].y[:]

        # compare horton - napmo
        print('spherical decomposition:', np.allclose(p_lm_n, p_lm_h))

        # Build rho again:

        # napmo
        rho_n = grid.atgrids[-1].evaluate_expansion(lmax, p_lm_n)
        print('rebuild napmo:', np.allclose(rho_n, rho))

        # horton
        rho_horton = molgrid.zeros()
        molgrid.eval_decomposition(
            density_decomposition, molgrid.centers[-1], rho_horton)
        print('rebuild horton:', np.allclose(rho_horton, rho))

        #############################
        # Integration test
        #############################

        # Solve Poisson
        ###############

        # horton
        print('horton')
        U_lm_h = solve_poisson_becke(density_decomposition)

        # napmo
        print('napmo')
        rho_n = p_ab(grid.points)
        U_lm_n = poisson_solver(grid, rho_n, lmax)

        # for i in range(lsize):
        #     plt.plot(r[::-1], U_lm_n[:, i], '-x', label='napmo'+str(i))
        #     plt.plot(rad, U_lm_h[i].y, '-x', label='horton'+str(i))
        # plt.legend()
        # plt.show()

        # Build the potential and calculate integral
        # ##########################################

        # napmo
        print('napmo')
        for i in range(grid.atgrids[-1].radial_grid.size):
            U_lm_n[i, :] /= grid.atgrids[-1].radial_grid.points[i]

        V_ab_n = grid.atgrids[-1].evaluate_expansion(lmax, U_lm_n)
        val_n = grid.integrate(rho_n * V_ab_n)
        n_ints.append(val_n)

        # horton
        print('horton')
        V_ab_h = molgrid.zeros()
        molgrid.eval_decomposition(U_lm_h, molgrid.centers[-1], V_ab_h)
        val_h = atgrid.integrate(rho, V_ab_h)
        h_ints.append(val_h)

        # print(V_ab_n)
        # plt.plot(r, V_ab_n[::angularPoints], '-x', label='napmo')
        # plt.plot(rad, V_ab_h[::angularPoints], '-x', label='horton')
        # plt.legend()
        # plt.show()

        # print('Potential:', np.allclose(V_ab_h, V_ab_n))
        print(val_h, val_n)

    plt.yscale('log')
    plt.xlabel('Radial Points')
    plt.ylabel('Error (log)')
    plt.plot(radialList, np.abs(np.array(h_ints) -
                                0.79788456080286518), label='horton')
    plt.plot(radialList, np.abs(np.array(n_ints) -
                                0.79788456080286518), label='napmo-7')

    plt.legend()
    plt.show()
    # plt.savefig('error.png')
