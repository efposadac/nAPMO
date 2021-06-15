# file: test_psi_numeric_becke_H2.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np

from ctypes import *


def test_psi_numeric_becke_H2(basis="STO-3G"):
    """
    test_psi_numeric_becke_H2
    * Closed shell
    * GTO: -1.117506 (CCCBDB, R=0.7122)
    * Limit: -1.133630
    * GRID: -1.1330638522586503 (100x110)
    * BASIS: STO-3G
    """
    file = ("""
  # Molecule definition:

  molecule {
    H 0.000 0.000  0.356
    H 0.000 0.000 -0.356
  }

  basis {
   e- """ + basis + """
  }

  scf {
    maxiter = 100
    hf
    numeric
  }

  grid {
    e- [100, 110]
  }
  """)

    # Build the system
    data = napmo.InputParser(file)
    system = napmo.NAPMO(data, pprint=False)

    # Build the solver
    solver = napmo.HF(system.system, options=system.data.scf, pprint=False)

    # Simulate SCF iteration (_a GTO-based, _n Grid-based)
    print("\n*** GTO-Based initialization")
    psi_a = solver.PSI[-1]
    ene_a = solver.scf.compute_energy_single(psi_a, show=True)

    # Compare to CCCBDB (GTO-based result)
    assert np.allclose(ene_a, -1.117506, rtol=1e-5)
    print("Energy GTO : {0:>18.14f}".format(ene_a))

    print("\n*** Grid-Based initialization")
    psi_n = solver.NPSI[-1]

    # psi_a.plot_dens(grid=psi_n.grid, kind="GTO-BASED-" + basis, xlim=[-1.0, 1.0], marker='-')

    # Calculate 1-body matrix
    psi_n.initialize()
    psi_n.compute_1body()

    # Calculate 2-body matrix
    psi_n.compute_2body(solver.get('direct'))

    # build Fock and compare
    psi_n.build_fock()

    # Iterate Once
    napmo.cext.wavefunction_iterate(byref(psi_n))

    # Energy
    energy = solver.scf.compute_energy_single(psi_n, show=True)
    print("Energy num : {0:>18.14f}".format(energy))

    Fgrid = analyse_psi(psi_n)

    # Aux basis
    E_tot = psi_optimization(psi_n, psi_a)
    # assert np.allclose(E_tot, -1.133630, rtol=1e-4)

    # psi_n.plot_dens(psi=psi_n.psi, kind="GRID-BASED-" + basis, xlim=[-1.0, 1.0], marker='--')

    return 0.0  # E_tot


def psi_optimization(psi, psi_a):
    print("\n*** Becke Opt Test")

    print("\n--- Convergence Test", psi._eta)

    converged = False
    E_prev = 0.0
    it = 1
    psi._debug = False
    prev_orb = psi_a.O
    while not converged:
        psi.optimize_psi()
        psi.compute_2body()
        psi.build_fock()

        napmo.cext.wavefunction_iterate(byref(psi))

        E_tot = psi.energy + psi.pce
        E_diff = np.abs(E_tot - E_prev)

        print("{0:>4d} {1:>18.14f} {2:>18.14f} {3:>18.14f} {4:>18.14f}".format(
            it, E_tot, E_diff, psi._optimize.delta_e.sum(), psi._optimize.delta_orb.sum()))

        if E_diff < 1.0e-7:
            converged = True

        E_prev = E_tot
        it += 1

    return E_tot


def analyse_psi(psi):
    # Fock operator
    T_num = (np.array([psi._grid.integrate(p * t)
                       for p, t in zip(psi.psi, psi.Tgrid)]) * psi.D.T).sum()

    V_num = (np.array([psi._grid.integrate(p * v)
                       for p, v in zip(psi.psi, psi.Vgrid)]) * psi.D.T).sum()

    aux_j = np.array([psi._grid.integrate(p * p * psi.Jpot) for p in psi.psi])

    J_num = (aux_j * 0.5 * psi.D.T).sum()

    aux_k = np.array([psi._grid.integrate(p * k)
                      for p, k in zip(psi.psi, psi.Kpot)])

    K_num = (psi._x_factor * aux_k * 0.5 * psi.D.T).sum()

    G_num = J_num + K_num

    # Fock representation in the grid (KET) ===  F = T + V + J - K
    Fgrid = np.array(
        [t + (v) + (
            0.5 * ((p * psi.Jpot) + (k * psi._x_factor)))
            for t, v, p, k in zip(psi.Tgrid, psi.Vgrid, psi.psi, psi.Kpot)])

    # Calculate (BRA)
    F_n = psi.get_operator_matrix(Fgrid)

    E_num = (psi.D.T * F_n).sum() + psi._pce

    # Test Orbitals
    aux_e = np.array([e * p for e, p in zip(psi.O, psi.psi)])
    E_n = psi.get_operator_matrix(aux_e)
    E_o_A = (E_n * psi.D.T).sum() + psi._pce
    E_o_B = psi.D.dot(psi.O).sum()
    E_t_o_A = 0.5 * (T_num + V_num + E_o_A)
    E_t_o_B = 0.5 * (T_num + V_num + E_o_B)

    print("\n--- PSI ANALYSIS")
    print("T num      : {0:>18.14f}".format(T_num))
    print("V num      : {0:>18.14f}".format(V_num))
    print("J num      : {0:>18.14f}".format(J_num))
    print("K num      : {0:>18.14f}".format(K_num))
    print("2body num  : {0:>18.14f}".format(G_num))
    print("Orbitals A : {0:>18.14f}".format(E_o_A))
    print("Orbitals B : {0:>18.14f}".format(E_o_B))
    print("Orbitals EA: {0:>18.14f}".format(E_t_o_A))
    print("Orbitals EB: {0:>18.14f}".format(E_t_o_B))
    print("Energy num : {0:>18.14f}".format(E_num))
    print("Energy SCF : {0:>18.14f}".format(psi._energy + psi._pce))

    return Fgrid


if __name__ == '__main__':
    test_psi_numeric_becke_H2()
    # energies = []
    # basis_list = ['STO-3G', 'CC-PVDZ', 'CC-PVTZ', 'CC-PVQZ', 'CC-PV5Z']
    # for basis in basis_list:
    #     energies.append(test_psi_numeric_becke_H2(basis))

    # for basis, energy in zip(basis_list, energies):
    #     print("Basis: %10s  Energy: %20.15f" % (basis, energy))

    # plt.show()
