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
        self.PSI = [napmo.PSIA(self.system.get_species(i),
                               self.system.point_charges)
                    for i in range(self.system.size_species)]

        self._pce = sum([psi.pce for psi in self.PSI])

        # SCF solver
        self.scf = napmo.SCF(options=self.options, pce=self._pce)

        # Numerical initialization
        self._mgrid = []
        self.NPSI = []

        if self.get('kind') is 'numeric':

            # Perform a single iteration (needed for numerical initialization)
            if len(self.PSI) > 1:
                self.scf.multi(self.PSI, pprint=True)
            else:
                self.scf.single(self.PSI[-1], pprint=True)

            # Build grids and numerical wavefunctions
            for p, key in enumerate(system.keys()):

                particle = system.get(key, {})
                if particle.get('is_electron'):
                    key = 'e-'

                aux = self.get('grid').get(key, None)

                if aux is None:
                    napmo.raise_exception(
                        ValueError,
                        "Grid specification not found!",
                        'Check the "grid" block in your input file ' + key + ' has not grid specifications')

                self._mgrid.append(napmo.BeckeGrid(particle,
                                                   aux.get('nrad', 100),
                                                   aux.get('nang', 110),
                                                   rtransform=aux.get('rtransform', None)))

                self._mgrid[-1].show()

                if self.get('hybrid', {}).get(self.PSI[p].symbol, 'N') == 'N':
                    self.NPSI.append(napmo.PSIN(self.PSI[p], self._mgrid[-1]))

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

        self._energy = self.scf._energy

        self.scf.show_results(self.PSI)

    def compute_numeric(self, pprint=True):
        """
        Numerical HF
        """

        # Multi-species
        if len(self.PSI) > 1:
            with napmo.runtime.timeblock('SCF Multi'):
                self.scf.nmulti(
                    self.NPSI, pprint=pprint)

        # Single-species
        else:
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
