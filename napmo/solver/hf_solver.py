# file: hf_solver.py
# nAPMO package
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

from ctypes import *
import numpy as np
import napmo


class HF(object):
    """
    Hartree-Fock solver
    """

    def __init__(self, system, options=None, pprint=True):
        super(HF, self).__init__()

        self.system = system

        self.options = {'method': 'hf',
                        'kind': 'analytic',
                        'grid': [100, 110]}

        if options:
            self.options.update(options)

        # Analytic initialization
        self.PSI = [napmo.PSIA(self.system.get_species(i),
                               self.system.point_charges,
                               self.system.total_mass,
                               options=self.options)
                    for i in range(self.system.size_species)]

        if 'hybrid' in self.options:
            self.options['kind'] = 'numeric'
            self.options['hybrid'] = {k: v for d in self.options.get('hybrid', {})
                                      for k, v in d.items()}

        self._pce = sum([psi.pce for psi in self.PSI])

        # SCF solver
        self.scf = napmo.SCF(options=self.options,
                             pce=self._pce, pprint=pprint)

        # Numerical initialization
        self._mgrid = []
        self.NPSI = []
        self.HPSI = []

        if self.get('kind') == 'numeric':

            # Perform a single iteration (needed for numerical initialization)
            if len(self.PSI) > 1:
                self.scf.multi(self.PSI, pprint=pprint)
            else:
                self.scf.single(self.PSI[-1], pprint=pprint)

            # Build grids and numerical wave-functions
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

                if not particle.get('is_nuclei', False):
                    self._mgrid.append(
                        napmo.BeckeGrid(
                            particle,
                            aux.get('nrad', 100),
                            aux.get('nang', 110),
                            rtransform=aux.get('rtransform', None),
                            file=aux.get('file', None),
                            ablmax=self.get('aux_basis_lmax', None),
                            abldep=self.get('aux_basis_ldep', None)
                        )
                    )

                    if particle.get('is_electron'):
                        electronic_grid = self._mgrid[-1]
                else:
                    if isinstance(electronic_grid, napmo.BeckeGrid):
                        self._mgrid.append(electronic_grid)
                        self._mgrid[-1]._symbol = particle.get('symbol')
                    else:
                        napmo.raise_exception(
                            ValueError,
                            "Grid specification not found!",
                            'Check the "Geometry" block in your input file and define electrons before nuclei')

                if pprint:
                    self._mgrid[-1].show()

                if 'hybrid' in self.options:

                    if self.get('hybrid', {}).get(self.PSI[p].symbol, 'N') == 'N':
                        self.HPSI.append(napmo.PSIN(
                            self.PSI[p], self._mgrid[0], debug=self.get('debug')))

                    elif self.get('hybrid', {}).get(self.PSI[p].symbol, 'N') == 'A':
                        self.HPSI.append(napmo.PSIH(self.system.get_species(p),
                                                    self.system.point_charges,
                                                    self.system.total_mass,
                                                    self._mgrid[0],
                                                    self.PSI[p]))

                else:
                    if self.get('hybrid', {}).get(self.PSI[p].symbol, 'N') == 'N':
                        self.NPSI.append(napmo.PSIN(
                            self.PSI[p], self._mgrid[-1],
                            debug=self.get('debug'),
                            aux_basis=self.get('aux_basis', False)
                        ))

        self._energy = 0.0

    def compute(self, pprint=True):
        """
        Computes APMO-HF or HF depending on the number of quantum species

        Args:
            pprint (bool): Whether to print or not the progress of the calculation.
        """
        if self.get('kind') == 'numeric':
            if 'hybrid' in self.options:
                self.compute_hybrid(pprint)
            else:
                self.compute_numeric(pprint)
        else:
            self.compute_analytic(pprint)

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

        if pprint:
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

        if pprint:
            self.scf.show_results(self.NPSI)

    def compute_hybrid(self, pprint=True):
        """
        Hybrid HF
        """

        # Only multi-species is expected.
        with napmo.runtime.timeblock('SCF Multi'):
            self.scf.hmulti(self.HPSI, pprint=pprint)

        self._energy = self.scf.energy

        if pprint:
            self.scf.show_results(self.HPSI)

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
