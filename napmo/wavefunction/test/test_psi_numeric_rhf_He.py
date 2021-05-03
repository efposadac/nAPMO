# file: test_psi_numeric.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

import napmo
import numpy as np

from ctypes import *
# import matplotlib.pyplot as plt


def test_psi_numeric_rhf_He(basis='STO-3G'):
    """
    * Closed shell
    * Single center
    * Single orbital
    * GTO: -2.807784 (CCCBDB)
    * GRID: -2.861679541872 (50x110)
    * BASIS: STO-3G
    """
    file = ("""
  # Molecule definition:

  molecule [e-:0:1] {
    He    0.000000  0.00000  0.0000
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
    e- [50, 110]
  }
  """)

    # Build the system
    data = napmo.InputParser(file)
    system = napmo.NAPMO(data, pprint=False)

    # Build the solver
    solver = napmo.HF(system.system, options=system.data.scf, pprint=False)

    # Simulate SCF iteration (_a GTO-based, _n Grid-based)
    psi_a = solver.PSI[-1]

    print("\n*** GTO-Based initialization")
    ene_a = solver.scf.compute_energy_single(psi_a, show=True)

    psi_n = solver.NPSI[-1]
    assert psi_n.species.get('size') == 2  # 2 electrons

    # psi_a.plot_dens(grid=psi_n.grid, kind="GTO-BASED-" + basis, xlim=[-1.0, 1.0], marker='-')

    # Calculate 1-body matrix
    psi_n.initialize()
    psi_n.compute_1body()

    # Operators H, T, V (H-core)
    assert np.allclose(psi_a.T, psi_n.T)
    assert np.allclose(psi_a.V, psi_n.V)
    assert np.allclose(psi_a.H, psi_n.H)

    # Calculate 2-body matrix
    psi_n.compute_2body(solver.get('direct'))
    assert np.allclose(psi_a.G, psi_n.G, rtol=1e-4)

    # build Fock and compare
    psi_n.build_fock()
    assert np.allclose(psi_a.F, psi_n.F, rtol=1e-4)

    # Iterate Once and Compare
    napmo.cext.wavefunction_iterate(byref(psi_n))
    assert np.allclose(psi_a.O, psi_n.O, rtol=1e-4)

    # Energy
    energy = psi_n._energy + psi_n.pce
    assert np.allclose(energy, psi_a._energy, rtol=1e-4)

    print("\n*** Grid-Based initialization")
    Fgrid = analyse_psi(psi_n)

    # Aux basis
    E_tot = aux_basis(psi_n)
    assert np.allclose(E_tot, -2.861679541872, rtol=1e-5)

    # psi_n.plot_dens(psi=psi_n.psi, kind="GRID-BASED-" + basis, xlim=[-1.0, 1.0], marker='--')

    return E_tot


def aux_basis(psi):
    print("\n*** Aux Basis Test")

    print("\n--- Convergence Test", psi._eta)

    converged = False
    E_prev = 0.0
    it = 1
    psi._debug = False
    while not converged:
        psi.optimize_psi()
        psi.compute_2body()
        psi.build_fock()

        napmo.cext.wavefunction_iterate(byref(psi))

        E_tot = psi.energy

        print(it, E_tot)

        if np.abs(E_tot - E_prev) < 1.0e-8:
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
    print("Energy SCF : {0:>18.14f}".format(psi._energy))

    return Fgrid


if __name__ == '__main__':
    test_psi_numeric_rhf_He()
    # energies = []
    # basis_list = ['STO-3G', 'CC-PVDZ', 'CC-PVTZ', 'CC-PVQZ', 'CC-PV5Z', 'CC-PV6Z']
    # for basis in basis_list:
    #     energies.append(test_psi_numeric_rhf_He(basis))

    # for basis, energy in zip(basis_list, energies):
    #     print("Basis: %10s  Energy: %20.15f" % (basis, energy))

    # plt.show()
