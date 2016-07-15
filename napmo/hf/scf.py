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
        Computes APMO-HF or HF depending on the number of quantum species

        Args:
            pprint (bool): Whether to print or not the progress of the calculation.
        """
        self.compute_1body_ints()
        self.compute_hcore()
        self.compute_guess()

        if len(self.PSI) > 1:
            with napmo.runtime.timeblock('SCF Multi'):
                self.iterate_multi()
        else:
            with napmo.runtime.timeblock('SCF Single'):
                self.iterate_single(self.PSI[-1], pprint=pprint)

        self.show_results()

        return self._energy

    def compute_1body_ints(self):
        """
        Computes one boy molecular integrals for all quantum species in the system, including:

        * Overlap
        * Kinetic
        * Nuclear potential
        """
        # Calculate Integrals
        with napmo.runtime.timeblock('1 body ints'):
            for psi in self.PSI:
                psi.compute_overlap()
                psi.compute_kinetic()
                psi.compute_nuclear()

    def compute_hcore(self):
        """
        Builds the one particle hcore.
        """
        for psi in self.PSI:
            psi.compute_hcore()

    def compute_guess(self):
        """
        Computes the density guess using the *HCORE* method:
        :math:`H C = e S C`
        """
        # Compute guess Density (HCORE)
        for psi in self.PSI:
            psi.compute_guess()

    def iterate_single(self, psi, pprint=False):
        """
        Perform SCF procedure for one species.

        Args:
            psi (WaveFunction) : WaveFunction object for one species.
            pprint (bool): Whether to print or not the progress of the calculation.
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

            with napmo.runtime.timeblock('DIIS'):
                napmo.cext.LibintInterface_diis(psi._diis, byref(psi))

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
        Perform SCF iteration for all species in this object.
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

                with napmo.runtime.timeblock('DIIS'):
                    napmo.cext.LibintInterface_diis(psi._diis, byref(psi))

            self.compute_energy()

            e_diff = self._energy - e_last

            print('{0:<4d} {1:>12.7f} {2:>12.7f} {3:>12.7f}'.
                  format(iterations, self._energy - self.pce, self._energy, e_diff))

    def compute_energy(self):
        """
        Computes the total energy for a multi-species system.
        """

        self._energy = 0.0
        self._coupling_energy = 0.0

        for psi in self.PSI:
            # Add independient particle energies
            self._energy += (psi.D * (psi.H + (0.5 * psi.G))).sum()

            # Calculate coupling energy
            self._coupling_energy += 0.5 * (psi.D * psi.J).sum()

        # Add point charges energy
        self._energy += self.pce

        # Add coupling Energy
        self._energy += self._coupling_energy

    def compute_energy_components(self):
        """
        Computes the total energy for a multi-species system.
        """
        self._kinetic_energy = 0.0
        self._coupling_energy = 0.0
        self._2body_energy = 0.0
        self._1body_energy = 0.0
        self._pcqui_energy = 0.0

        for psi in self.PSI:
            # Calculate kinetic energy
            self._kinetic_energy += (psi.D * psi.T).sum()

            # Calculate kinetic energy
            self._1body_energy += (psi.D * psi.H).sum()

            # Calculate point charges - quantum interaction energy
            self._pcqui_energy += (psi.D * psi.H).sum() - (psi.D * psi.T).sum()

            # Calculate repulsion energy
            self._2body_energy += 0.5 * (psi.D * psi.G).sum()

            # Calculate coupling energy
            self._coupling_energy += 0.5 * (psi.D * psi.J).sum()

        # Calculate potential energy
        self._potential_energy = (self.pce +
                                  self._2body_energy +
                                  self._pcqui_energy +
                                  self._coupling_energy)

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

    def show_results(self):
        """
        Prints detailed results of the calculation
        """
        self.compute_energy_components()

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
                          -----------------
  Potential energy         {8:>16.11f}

""".format(self._potential_energy,
           self._kinetic_energy,
           self._potential_energy + self._kinetic_energy,
           self._potential_energy / self._kinetic_energy,
           self.pce,
           self._pcqui_energy,
           self._2body_energy,
           self._coupling_energy,
           self._potential_energy))

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
