# file: test_psi_numeric_ab_Li.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np

from ctypes import *


def test_psi_numeric_ab_Li(basis="STO-3G"):
    """
    test_psi_numeric_ab_Li
    * Open shell
    * Single center
    * GTO: -7.315526 (CCCBDB)
    * HF CBSE: -7.432720216
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
    aux_basis
    aux_basis_ldep = 1.0e-6
    aux_basis_lmax = 0
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
    test_psi_numeric_ab_Li()
    # energies = []
    # basis_list = ['STO-3G', 'CC-PVDZ', 'CC-PVTZ', 'CC-PVQZ', 'CC-PV5Z']
    # for basis in basis_list:
    #     energies.append(test_psi_numeric_ab_Li(basis))

    # for basis, energy in zip(basis_list, energies):
    #     print("Basis: %10s  Energy: %20.15f" % (basis, energy))

    # plt.show()
