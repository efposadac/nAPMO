# file: test_psi_numeric_ab_PsH.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np

from ctypes import *


def test_psi_numeric_ab_PsH(basis="AUG-CC-PVTZ"):
    """
    test_psi_numeric_ab_PsH
    * APMO
    * GTO: -0.6662409
    * HF CBSE: -0.6670846504
    * GRID:  (100x110)
    * BASIS: aug-cc-pVTZ
    """
    file = ("""
  # Molecule definition:

  molecule [e-:1] {
    H  0.00 0.00 0.00
    E+ 0.00 0.00 0.00
  }

  basis {
   e- """ + basis + """
   E+ E+-H-7SP-AUG-CC-PVTZ
  }

  scf {
    maxiter = 100
    hf
    numeric
    aux_basis
    aux_basis_ldep = 1.0e-6
    aux_basis_lmax = 0
  }

  grid {
    e- [0.5858151015155881, 100, 50]
    E+ [0.5858151015155881, 100, 50]
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

    print("Energy GTO : {0:>18.14f}".format(ene_a))

    # Compare to GTO-Based (GTO-based result)
    assert np.allclose(ene_a, -0.6662409, rtol=1e-5)

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
    E_tot = psi_optimization(psi_n)
    assert np.allclose(E_tot, -0.6670846504, rtol=1e-3)

    # psi_n.plot_dens(psi=psi_n.psi, kind="GRID-BASED-" + basis, xlim=[-1.0, 1.0], marker='--')

    return E_tot


def psi_optimization(psi):
    print("\n*** Aux Basis Opt Test")
    converged = False
    E_prev = 0.0
    it = 1
    # for i in range(4):
    while not converged:
        for p in psi:
            if it > 1:
                p.optimize_psi(other_psi=psi)

            p.compute_2body()
            p.build_fock()
            napmo.cext.wavefunction_iterate(byref(p))

        for p in psi:
            p.compute_coupling(psi)

        E_tot = compute_energy(psi)
        E_diff = E_tot - E_prev

        print(it, E_tot, E_diff)

        if np.abs(E_diff) < 1.0e-7:
            converged = True

        E_prev = E_tot
        it += 1

    return E_tot


def compute_energy(PSI):
    """
    Computes the total energy for a multi-species system.
    """

    energy = 0.0
    coupling_energy = 0.0
    xc_energy = 0.0

    for psi in PSI:
        # Add independent particle energies
        aux1 = (psi.D.T * (psi.H + (0.5 * psi.G))).sum()
        energy += aux1

        # Calculate coupling energy
        aux2 = 0.5 * (psi.D.T * psi.J).sum()
        coupling_energy += aux2

        # Calculate XC energy
        xc_energy += psi._xc_energy

        # print("\n--- ENERGY ANALYSIS", psi.symbol)
        # print("Single     : {0:>18.14f}".format(aux1))
        # print("Coupling   : {0:>18.14f}".format(aux2))
        # print("xc         : {0:>18.14f}".format(psi._xc_energy))

    # Add point charges energy
    energy += psi._pce

    # Add coupling Energy
    energy += coupling_energy

    # Add XC energy
    energy += xc_energy

    return energy


if __name__ == '__main__':
    test_psi_numeric_ab_PsH()
    # energies = []
    # basis_list = ['STO-3G', 'CC-PVDZ', 'CC-PVTZ', 'CC-PVQZ', 'CC-PV5Z']
    # for basis in basis_list:
    #     energies.append(test_psi_numeric_ab_PsH(basis))

    # for basis, energy in zip(basis_list, energies):
    #     print("Basis: %10s  Energy: %20.15f" % (basis, energy))

    # plt.show()
