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


class SCF(object):
    """
    Hartree-Fock Implementation
    """

    def __init__(self, system, options=None):
        super(SCF, self).__init__()
        print('\nStarting SCF Calculation...')

        self.options = {'maxiter': 100, 'eps_e': 1e-7,  'eps_d': 1e-7,
                        'method': 'hf', 'kind': 'analytic', 'direct': False}

        if options:
            self.options.update(options)

        print(self)

        self.system = system

        # Initialize WaveFunction:
        self.PSI = [
            napmo.WaveFunction(self.system.get_species(i),
                               self.system.point_charges)
            for i in range(self.system.size_species)]

        self._energy = 0.0
        self._pce = sum([psi.pce for psi in self.PSI])

        print("Point charges energy: {0:<12.8f}".format(self._pce))

    def compute(self, pprint=True):
        """
        """
        self.compute_1body_ints()
        self.compute_hamiltonian()
        self.compute_guess()

        if len(self.PSI) > 1:
            with napmo.runtime.timeblock('SCF Multi'):
                self.iterate_multi()
        else:
            with napmo.runtime.timeblock('SCF Single'):
                self.iterate_single(self.PSI[-1], pprint=pprint)

        return self._energy

    def compute_1body_ints(self):
        """
        """
        # Calculate Integrals
        with napmo.runtime.timeblock('1 body ints'):
            for psi in self.PSI:
                psi.compute_overlap()
                psi.compute_kinetic()
                psi.compute_nuclear()

    def compute_hamiltonian(self):
        """
        """
        for psi in self.PSI:
            psi.compute_hamiltonian()

    def compute_guess(self):
        """
        """
        # Compute guess Density (HCORE)
        for psi in self.PSI:
            psi.compute_guess()

    def iterate_single(self, psi, pprint=False):
        """
        """

        if pprint:
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

            psi.build_fock()

            # solve F C = e S C
            with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                napmo.cext.wavefunction_iterate(byref(psi))

            e_diff = psi._energy - e_last

            self._energy = psi._energy + psi.pce

            # print results
            if pprint:
                print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f} {4:>12.7f}'.
                      format(iterations, psi._energy, self._energy, e_diff, psi._rmsd))

    def iterate_multi(self):
        """
        """

        print('{0:5s}  {1:^10s} {2:>12s} {3:>12s}'
              .format("\nIter", "Energy", "Total E", "Delta(E)"))

        iterations = 0
        e_diff = 1

        while (iterations < self.get('maxiter') and
               np.abs(e_diff) > self.get('eps_e')):

            iterations += 1

            e_last = self._energy

            for psi in self.PSI:
                # Calculate 2 body Matrix
                with napmo.runtime.timeblock('2 body ints'):
                    psi.compute_2body(self.get('direct'))

                psi.build_fock()

                # solve F C = e S C
                with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                    napmo.cext.wavefunction_iterate(byref(psi))

                with napmo.runtime.timeblock('Coupling ints'):
                    psi.compute_couling(self.PSI, direct=self.get('direct'))

                psi.build_fock()

                # solve F C = e S C
                with napmo.runtime.timeblock('Self-Adjoint eigen solver'):
                    napmo.cext.wavefunction_iterate(byref(psi))

            self.compute_energy()

            e_diff = self._energy - e_last

            print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f}'.
                  format(iterations, self._energy - self.pce, self._energy, e_diff))

    def compute_energy(self):
        """
        """

        self._energy = 0.0
        self._coupling_energy = 0.0

        for psi in self.PSI:
            # Add independient particle energies
            self._energy += psi._energy

            # Calculate coupling energy
            self._coupling_energy += 0.5 * (psi.D * psi.J).sum()

        # print("Energy species: ", self._energy)
        # print("Coupling Energy: ", self._coupling_energy)
        # print("Point charges: ", self.pce)

        # Add point charges energy
        self._energy += self.pce

        # Add coupling Energy
        self._energy += self._coupling_energy

        # print("Total: ", self._energy)
        # Add electronic repulsion energy (open-shell case)
        # Not implemented

    def get(self, key, default=None):
        """
        """
        return self.options.get(key, default)

    @property
    def pce(self):
        return self._pce

    def __repr__(self):
        out = ("""\nSCF setup:

Method:   {0:<10s}
Kind:     {1:<10s}
Direct:   {2:<10s}
E Tol:    {3:<10.3e}
Dens Tol: {4:<10.3e}
""".format(self.get('method'),
           self.get('kind'),
           str(self.get('direct')),
           self.get('eps_e'),
           self.get('eps_d')
           ))

        return out
