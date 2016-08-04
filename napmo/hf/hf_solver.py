from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


class HF(object):
    """
    Hartree-Fock solver
    """

    def __init__(self, system, options=None):
        super(HF, self).__init__()

        self.system = system

        self.options = {'method': 'hf',
                        'kind': 'analytic',
                        'grid': [50, 110]}

        if options:
            self.options.update(options)

        if 'hybrid' in self.options:
            self.options['kind'] = 'numeric'
            self.options['hybrid'] = {k: v for d in self.options.get('hybrid', {})
                                      for k, v in d.items()}

        # Analytic initialization
        self.PSI = [napmo.WaveFunction(self.system.get_species(i),
                                       self.system.point_charges)
                    for i in range(self.system.size_species)]

        self._pce = sum([psi.pce for psi in self.PSI])

        # SCF solver
        self.scf = napmo.SCF(options=self.options, pce=self._pce)

        # Perform a single iteration (needed for numerical initialization)
        for psi in self.PSI:
            self.scf.iteration(psi)

        # Numerical initialization
        if self.get('kind') is 'numeric':

            # Grid definition (Only atoms for now)
            self._mgrid = napmo.BeckeGrid(system.get(
                'e-'), self.get('grid')[0], self.get('grid')[1])

            self._mgrid.show()

            self.NPSI = [napmo.NWaveFunction(psi, self._mgrid)
                         for psi in self.PSI
                         if self.get('hybrid', {}).get(psi.symbol, 'N') == 'N']

        self._energy = 0.0

    def compute(self, pprint=True):
        """
        Computes APMO-HF or HF depending on the number of quantum species

        Args:
            pprint (bool): Whether to print or not the progress of the calculation.
        """
        if self.get('kind') is 'numeric':
            if 'hybrid' in self.options:
                self.compute_hybrid()
            else:
                self.compute_numeric()
        else:
            self.compute_analytic()

        return self._energy

    def compute_analytic(self, pprint=True):
        """
        Conventional HF using GTO's
        """

        # Multi-species
        if len(self.PSI) > 1:
            with napmo.runtime.timeblock('SCF Multi'):
                self.scf.multi(
                    self.PSI, pprint=pprint)

        # Single-species
        else:
            with napmo.runtime.timeblock('SCF Single'):
                self.scf.single(
                    self.PSI[-1], pprint=pprint)

        self._energy = self.scf.energy

        self.scf.show_results(self.PSI)

    def compute_numeric(self, pprint=True):
        with napmo.runtime.timeblock('SCF Numeric'):
            self.scf.nsingle(
                self.NPSI[-1], pprint=pprint)

        self._energy = self.scf.energy

        self.scf.show_results(self.NPSI)

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
