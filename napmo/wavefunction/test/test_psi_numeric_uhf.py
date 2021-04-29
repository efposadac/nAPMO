# file: test_psi_numeric_uhf.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

import napmo
import numpy as np

from ctypes import *


def test_psi_numeric_uhf(basis="STO-3G"):
    """
    * Open shell
    * Single center
    * GTO: -7.315526 (CCCBDB)
    * HF CBSE: -7.432720586
    * GRID: -7.4327362568583855 (50x14)
    * BASIS: STO-3G
    """
    file = ("""
  # Molecule definition:

  molecule {
    Li    0.000000  0.00000  0.0000
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
    e- [50, 14]
  }
  """)

    # Build the system
    data = napmo.InputParser(file)
    system = napmo.NAPMO(data, pprint=True)

    # Build the solver
    solver = napmo.HF(system.system, options=system.data.scf, pprint=True)

    # Simulate SCF iteration (_a GTO-based, _n Grid-based)
    print("\n*** GTO-Based initialization")
    psi_a = solver.PSI
    solver.scf.multi(psi_a, pprint=False)
    ene_a = solver.scf._energy

    # Compare to CCCBDB (GTO-based result)
    assert np.allclose(ene_a, -7.315526, rtol=1e-5)

    print("\n*** Grid-Based initialization")
    psi_n = solver.NPSI

    # psi_a.plot_dens(grid=psi_n.grid, kind="GTO-BASED-" + basis, xlim=[-1.0, 1.0], marker='-')

    for psi in psi_n:
        # Calculate 1-body matrix
        psi.initialize()
        psi.compute_1body()

        # Calculate 2-body matrix
        psi.compute_2body(solver.get('direct'))

        # build Fock and compare
        psi.build_fock()

        # Iterate Once
        napmo.cext.wavefunction_iterate(byref(psi))

    for psi in psi_n:
        # Add coupling
        psi.compute_coupling(psi_n)

    # Energy
    energy = compute_energy(psi_n)
    assert np.allclose(energy, ene_a, rtol=1e-4)
    print("Energy num : {0:>18.14f}".format(energy))

    # Aux basis
    E_tot = aux_basis(psi_n)
    assert np.allclose(E_tot, -7.4327362568583855, rtol=1e-4)

    # psi_n.plot_dens(psi=psi_n.psi, kind="GRID-BASED-" + basis, xlim=[-1.0, 1.0], marker='--')

    return E_tot


def aux_basis(psi):
    print("\n*** Aux Basis Test")

    print("\n--- Convergence Test", psi[-1]._eta)

    converged = False
    E_prev = 0.0
    it = 1
    # psi._debug = False
    while not converged:
        for p in psi:
            if it > 1:
                p.optimize_psi(other_psi=psi)
            p.compute_2body()
            p.build_fock()

            napmo.cext.wavefunction_iterate(byref(p))
            p.compute_coupling(psi)

            E_tot = p.energy

        E_tot = compute_energy(psi)
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


def compute_energy(PSI):
    """
    Computes the total energy for a multi-species system.
    """

    energy = 0.0
    coupling_energy = 0.0
    xc_energy = 0.0

    for psi in PSI:
        # Add independent particle energies
        energy += (psi.D.T * (psi.H + (0.5 * psi.G))).sum()

        # Calculate coupling energy
        coupling_energy += 0.5 * (psi.D.T * psi.J).sum()

        # Calculate XC energy
        xc_energy += psi._xc_energy

    # Add point charges energy
    energy += psi._pce

    # Add coupling Energy
    energy += coupling_energy

    # Add XC energy
    energy += xc_energy

    return energy


if __name__ == '__main__':
    test_psi_numeric_uhf()
    # energies = []
    # basis_list = ['STO-3G', 'CC-PVDZ', 'CC-PVTZ', 'CC-PVQZ', 'CC-PV5Z']
    # for basis in basis_list:
    #     energies.append(test_psi_numeric_uhf(basis))

    # for basis, energy in zip(basis_list, energies):
    #     print("Basis: %10s  Energy: %20.15f" % (basis, energy))

    # plt.show()
