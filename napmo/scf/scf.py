# file: scf.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo
import matplotlib.pyplot as plt

import sys
from scipy import linalg as SLA


class SCF(object):
    """
    SCF Implementation

    Args:
        options (dict): Options to handle the SCF calculation
        pce (double) : Point charges energy
    """

    def __init__(self, options=None, pce=0.0, pprint=True):
        super(SCF, self).__init__()
        self.options = {'maxiter': 100,
                        'eps_e': 1e-7,
                        'eps_n': 1e-6,
                        'eps_d': 1e-6,
                        'eps_r': 1e-4,
                        'method': 'hf',
                        'kind': 'analytic',
                        'direct': False,
                        'print': True,
                        'debug': False}

        if options:
            self.options.update(options)

        if self.get('kind') is 'numeric':
            self.options['print'] = False

        if pprint:
            print(self)

        self._energy = 0.0
        self._pce = pce

        if pprint:
            print("Point charges energy: {0:<12.8f}".format(self._pce))

    def iteration_single(self, psi):
        """
        Performs a single iteration for one species.

        Args:
            psi (WaveFunction) : WaveFunction object for one species.
        """
        with napmo.runtime.timeblock('2 body ints'):
            psi.compute_2body(self.get('direct'))
            # psi.compute_exccor()

        psi.build_fock()

        # solve F C = e S C
        with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
            napmo.cext.wavefunction_iterate(byref(psi))

        self._energy = psi._energy + psi.pce

        if self.get('debug'):
            print('Single iteration {0:<3s} {1:>12.7f} {2:>12.7f} {3:>12.7f}'.
                  format(psi.symbol, psi._energy, self._energy, psi._rmsd))

    def single(self, psi, pprint=False, diis=True, other_psi=None):
        """
        Perform SCF procedure for one species.

        Args:
            psi (WaveFunction) : WaveFunction object for one species.
            pprint (bool): Whether to print or not the progress of the calculation.
        """

        if pprint:
            print('\nStarting Single SCF Calculation...')
            print('{0:5s}  {1:^10s} {2:>12s} {3:>12s} {4:>12s}'
                  .format("\nIter", "E (" + psi.symbol + ")", "Total E", "Delta(E)", "RMS(D)"))

        iterations = 0
        e_diff = 1

        while (iterations < self.get('maxiter') and
               np.abs(e_diff) > self.get('eps_e') and
               np.abs(psi._rmsd) > self.get('eps_d')):

            iterations += 1
            e_last = psi._energy

            with napmo.runtime.timeblock('2 body ints'):
                psi.compute_2body(self.get('direct'))
                psi.compute_exccor()

            if other_psi is not None:
                with napmo.runtime.timeblock('Coupling ints'):
                    psi.compute_coupling(other_psi, direct=self.get('direct'))
                    # psi.compute_cor2species(other_psi)

            psi.build_fock()

            # Iterate
            if iterations > 2:
                psi.F[:] = psi.convergence.damping(psi.F, psi.D)
            else:
                psi.convergence = napmo.Convergence(psi.F, psi.D)

            # if diis:
            #     with napmo.runtime.timeblock('DIIS'):
            #         napmo.cext.LibintInterface_diis(psi._diis, byref(psi))

            # solve F C = e S C
            with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                napmo.cext.wavefunction_iterate(byref(psi))

            e_diff = psi._energy - e_last

            self._energy = psi._energy + psi.pce

            # print results
            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f} {4:>12.7f}'.
                      format(iterations, psi._energy, self._energy, e_diff, psi._rmsd))

        if self.get('debug') and not isinstance(psi, napmo.PSIO):

            if not isinstance(psi, napmo.PSIO):
                grid = napmo.BeckeGrid(psi.species, 500, 110)
                psi.plot_dens(grid, kind="anal")
                plt.show()

    def multi(self, PSI, pprint=True, case=0):
        """
        Perform SCF iteration for all species in this object.

        Args:
            psi (WaveFunction) : WaveFunction object for one species.
            pprint (bool): Whether to print or not the progress of the calculation.
        """

        if pprint and case is 0:
            print('\nStarting Multi SCF Calculation...')
            print('{0:5s}  {1:^10s} {2:>12s} {3:>12s}'
                  .format("\nIter", "Energy", "Total E", "Delta(E)"))

        if pprint and case > 0:
            print('\nRestarting Multi SCF Calculation...')
            print('{0:5s}  {1:^10s} {2:>12s} {3:>12s}'
                  .format("\nIter", "Energy", "Total E", "Delta(E)"))

        iterations = 0
        e_diff = 1

        # Initialization
        for psi in PSI:
            # Calculate 2 body Matrix
            psi.compute_2body(self.get('direct'))
            # psi.compute_exccor()

        while (iterations < self.get('maxiter') and
               np.abs(e_diff) > self.get('eps_e')):

            iterations += 1

            e_last = self._energy

            if case is 0:

                for psi in PSI:
                    psi.build_fock()

                    # Iterate
                    if iterations > 2:
                        psi.F[:] = psi.convergence.damping(psi.F, psi.D)
                    else:
                        psi.convergence = napmo.Convergence(psi.F, psi.D)

                    with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                        napmo.cext.wavefunction_iterate(byref(psi))

                for psi in PSI:

                    # Calculate 2 body Matrix
                    with napmo.runtime.timeblock('2 body ints'):
                        psi.compute_2body(self.get('direct'))
                        # psi.compute_exccor()

                    with napmo.runtime.timeblock('Coupling ints'):
                        psi.compute_coupling(PSI, direct=self.get('direct'))
                        # psi.compute_cor2species(PSI)

            if case is 1:
                for psi in PSI:

                    with napmo.runtime.timeblock('2 body ints'):
                        psi.compute_2body(self.get('direct'))
                        # psi.compute_exccor()

                    psi.build_fock()

                    # solve F C = e S C
                    with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                        napmo.cext.wavefunction_iterate(byref(psi))

                    with napmo.runtime.timeblock('Coupling ints'):
                        psi.compute_coupling(PSI, direct=self.get('direct'))
                        # psi.compute_cor2species(PSI)

                    self.single(psi, pprint=False)

                    if self.get('debug'):
                        print("Single particle energy for " +
                              psi.symbol + ":", psi._energy + psi.pce)

            if case is 2:
                print('\n Calculation does not converge!')
                return

            self.compute_energy(PSI)

            e_diff = self._energy - e_last

            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f}'.
                      format(iterations, self._energy - self.pce, self._energy, e_diff))

        if iterations >= self.get('maxiter'):
            case += 1
            self.multi(PSI, case=case)
            return

        print('{0:11s} {1:>12.7f}'.format("\nANALYTICAL ", self._energy))

        if self.get('debug'):
            for psi in PSI:
                grid = napmo.BeckeGrid(psi.species, 100, 110)
                psi.plot_dens(grid, kind="anal")
                # plt.show()
                # plt.savefig('analytic_dens.png')

    def nsingle(self, psi, pprint=True):
        """
        Perform Numerical SCF procedure for one species.

        Args:
            psi (WaveFunction) : WaveFunction object for one species.
            pprint (bool): Whether to print or not the progress of the calculation.
        """

        if pprint:
            print('\nStarting Single NSCF Calculation...')
            print('{0:5s}  {1:^10s} {2:^12s} {3:^12s} {4:^12s}'
                  .format("\n***Iter", "E (" + psi.symbol + ")", "Total E", "Delta(E)", "Delta Orb"))

        iterations = 1
        e_diff = 1.0
        e_last = psi._energy

        converged = False
        while (not converged):

            if iterations > 1:
                # Compute \Psi (eq. 14) through conventional SCF for
                # \psi = a \phi + b \Delta \phi
                psi.optimize_psi(self)

            e_last = psi._energy

            with napmo.runtime.timeblock('Numerical 2 body'):
                psi.compute_2body(self.get('direct'))
                # psi.compute_exccor()

            psi.build_fock()

            # if iterations > 2:
            #     psi.F[:] = psi.convergence.damping(psi.F, psi.D)
            # else:
            #     psi.convergence = napmo.Convergence(psi.F, psi.D)

            # solve F C = e S C
            with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                napmo.cext.wavefunction_iterate(byref(psi))

            e_diff = psi._energy - e_last

            self._energy = psi._energy + psi.pce

            # print results
            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7g} {4:>12s}'.
                      format(iterations, psi._energy, self._energy, e_diff, str(psi._res)))

            if iterations > self.get('maxiter') or np.abs(e_diff) < self.get('eps_n'):
                if iterations > 1:
                    converged = True

            iterations += 1

        # if self.get('debug') and not isinstance(psi, napmo.PSIO):
        #     if not isinstance(psi, napmo.PSIO):
        #         psi.plot_dens(kind="num")
        #         plt.show()

    def nmulti(self, PSI, pprint=True):
        """
        Perform SCF iteration for all species in this object.
        """

        if pprint:
            print('\nStarting Multi NSCF Calculation...')
            print('{0:5s}  {1:^10s} {2:>12s} {3:>12s} {4:>12s}'
                  .format("\nIter", "Energy", "Total E", "Delta(E)", "Delta orb"))

        iterations = 1
        e_diff = 1
        converged = False
        while (not converged):

            e_last = self._energy

            if iterations > 1:
                for psi in PSI:
                    psi.optimize_psi(self, PSI)

            for psi in PSI:

                # Calculate 2 body Matrix
                with napmo.runtime.timeblock('2 body ints'):
                    psi.compute_2body(self.get('direct'))
                    # psi.compute_exccor()

                psi.build_fock()

                # if iterations > 2:
                #     psi.F[:] = psi.convergence.damping(psi.F, psi.D)
                # else:
                #     psi.convergence = napmo.Convergence(psi.F, psi.D)

                # solve F C = e S C
                with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                    napmo.cext.wavefunction_iterate(byref(psi))

                with napmo.runtime.timeblock('Coupling ints'):
                    psi.compute_coupling(PSI, direct=self.get('direct'))
                    # psi.compute_cor2species(PSI)

                if self.get('debug'):
                    print("Single particle energy for " +
                          psi.symbol + ":", psi._energy + psi.pce)

            self.compute_energy(PSI)

            e_diff = self._energy - e_last

            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7g} {4:>12s}'.
                      format(iterations, self._energy - self.pce, self._energy, e_diff,
                             str([str(psi._res) for psi in PSI])))

            if iterations > self.get('maxiter') or np.abs(e_diff) < self.get('eps_n'):
                if iterations > 1:
                    converged = True

            iterations += 1

        # if self.get('debug'):
        # for psi in PSI:
        #     psi.plot_dens(kind="num")
            # plt.show()
            # plt.savefig('numeric_dens.png')

    def hmulti(self, PSI, pprint=True):
        """
        Perform SCF iteration for all species in this object.
        """

        if pprint:
            print('\nStarting HYBRID Multi NSCF Calculation...')
            print('{0:5s}  {1:^10s} {2:>12s} {3:>12s}'
                  .format("\nIter", "Energy", "Total E", "Delta(E)"))

        iterations = 0
        e_diff = 1

        while (iterations < self.get('maxiter') and
               np.abs(e_diff) > self.get('eps_n')):

            # while (iterations < 3):
            iterations += 1

            e_last = self._energy

            if iterations > 1:
                for psi in PSI:
                    psi.optimize_psi(self, PSI)

            for psi in PSI:

                # Calculate 2 body Matrix
                psi.compute_2body(self.get('direct'))
                # psi.compute_exccor()
                psi.build_fock()

                if iterations > 2:
                    psi.F[:] = psi.convergence.damping(psi.F, psi.D)
                else:
                    psi.convergence = napmo.Convergence(psi.F, psi.D)

                # solve F C = e S C
                with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                    napmo.cext.wavefunction_iterate(byref(psi))

                with napmo.runtime.timeblock('Coupling ints'):
                    psi.compute_coupling(PSI, direct=self.get('direct'))
                    # psi.compute_cor2species(PSI)

                if psi.symbol == 'e-':
                    self.nsingle(psi, pprint=False)

                if self.get('debug'):
                    print(
                        "Single particle energy for " +
                        psi.symbol + ":", psi._energy + psi.pce)

            self.compute_energy(PSI)

            e_diff = self._energy - e_last

            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f} {4:>12s}'.
                      format(iterations, self._energy - self.pce, self._energy, e_diff, str([str(psi._rmsd) for psi in PSI])))

        if self.get('debug'):
            for psi in PSI:
                psi.plot_dens(kind="num")
                # plt.show()
                # plt.savefig('numeric_dens.png')

    def compute_energy(self, PSI):
        """
        Computes the total energy for a multi-species system.
        """

        self._energy = 0.0
        self._coupling_energy = 0.0

        for psi in PSI:
            # Add independient particle energies
            self._energy += (psi.D.T * (psi.H + (0.5 * psi.G))).sum()+psi._ecenergy

            # Calculate coupling energy
            self._coupling_energy += 0.5 * (psi.D.T * psi.J).sum()

            if self.get('debug'):
                print("coupling " + psi.symbol + ": ",
                      0.5 * (psi.D * psi.J).sum())

        # Add point charges energy
        self._energy += self.pce

        # Add coupling Energy
        self._energy += self._coupling_energy

    def compute_energy_components(self, PSI):
        """
        Computes the total energy for a multi-species system.
        """
        self._kinetic_energy = 0.0
        self._1body_energy = 0.0
        self._pcqui_energy = 0.0
        self._2body_energy = 0.0
        self._coupling_energy = 0.0
        self._ec_energy = 0.0

        for psi in PSI:
            # Calculate kinetic energy
            self._kinetic_energy += (psi.D * psi.T).sum()

            if self.get('debug'):
                print("Kinetic Energy " + psi.symbol +
                      ": ",  (psi.D * psi.T).sum())

            # Calculate one body energy
            self._1body_energy += (psi.D * psi.H).sum()

            # Calculate point charges - quantum interaction energy
            self._pcqui_energy += (psi.D * psi.V).sum()

            # Calculate repulsion energy
            self._2body_energy += 0.5 * (psi.D * psi.G).sum()

            # Calculate coupling energy
            self._coupling_energy += 0.5 * (psi.D * psi.J).sum()

            # Calculate exchange correlation energy
            self._ec_energy += psi._ecenergy

        # Calculate potential energy
        self._potential_energy = (self.pce +
                                  self._2body_energy +
                                  self._pcqui_energy +
                                  self._coupling_energy +
                                  self._ec_energy)

        self._energy = self._potential_energy + self._kinetic_energy

    def get(self, key, default=None):
        """
        Returns the option ``key`` of the SCF object
        """
        return self.options.get(key, default)

    @property
    def pce(self):
        """
        The point charges energy
        """
        return self._pce

    @property
    def energy(self):
        """
        The total energy
        """
        return self._energy

    def show_results(self, PSI):
        """
        Prints detailed results of the calculation
        """
        self.compute_energy_components(PSI)

        print("""\n
End SCF calculation
--------------------------------------------------

Hartree-Fock Results:
---------------------

  Total potential energy   {0:>16.11f}
  Total kinetic energy     {1:>16.11f}
                          -----------------
  Total Energy             {2:>16.11f}

  Virial ratio (V/T)       {3:>16.11f}

  Potential energy components:

  Point charges energy     {4:>16.11f}
  Quantum-Point energy     {5:>16.11f}
  Repulsion energy         {6:>16.11f}
  Coupling energy          {7:>16.11f}
  Exc. corr. energy        {8:>16.11f}
                          -----------------
  Potential energy         {9:>16.11f}

""".format(self._potential_energy,
           self._kinetic_energy,
           self._potential_energy + self._kinetic_energy,
           np.abs(self._potential_energy / self._kinetic_energy),
           self.pce,
           self._pcqui_energy,
           self._2body_energy,
           self._coupling_energy,
           self._ec_energy,
           self._potential_energy))

        for psi in PSI:
            print("Orbitals: ", psi.symbol, "\n", psi.O, "\n")
            print("Kinetic energy: ", psi.symbol, (psi.D * psi.T).sum())
            print("Quantum-Point energy: ", psi.symbol, (psi.D * psi.V).sum())
            print("Repulsion energy: ", psi.symbol, 0.5 * (psi.D * psi.G).sum())
            print("Coupling energy: ", psi.symbol, 0.5 * (psi.D * psi.J).sum())
            print("Exc. corr. energy (fix): ", psi.symbol, psi._ecenergy, "\n")

    def __repr__(self):
        out = ("""\nSCF setup:

Method:      {0:<10s}
Kind:        {1:<10s}
Direct:      {2:<10s}
E Tol:       {3:<10.3e}
Dens Tol:    {4:<10.3e}
E Tol NUM:   {5:<10.3e}
""".format(self.get('method'),
           self.get('kind'),
           str(self.get('direct')),
           self.get('eps_e'),
           self.get('eps_d'), 
           self.get('eps_n')
           ))

        return out