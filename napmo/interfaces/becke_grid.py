# file: becke_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from __future__ import division
from __future__ import print_function

import numpy as np
from ctypes import *

from napmo.interfaces.c_binding import *
from napmo.utilities.constants import *
from napmo.interfaces.atomic_grid import *


class BeckeGrid(Structure):
    """
    This class creates the Becke grid.

    References:
        Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
        n_radial (int, optional): Number of radial points. Default is 40
        n_angular (int, optional): Number of angular points. Default is 110
    """
    _fields_ = [
        ("n_radial", c_int),
        ("n_angular", c_int),
        ("_radial_abscissas", POINTER(c_double)),
        ("_radial_weights", POINTER(c_double)),
        ("_angular_theta", POINTER(c_double)),
        ("_angular_phi", POINTER(c_double)),
        ("_angular_weights", POINTER(c_double))
    ]

    def __init__(self, n_radial=40, n_angular=110):
        super(BeckeGrid, self).__init__()
        self.n_radial = c_int(n_radial)
        self.n_angular = c_int(n_angular)
        napmo_library.grid_init(byref(self))
        self.size = self.n_radial * self.n_angular
        self.expanded = False
        self.spherical = True
        self.xyz = np.empty((self.size, 3))
        self.w = np.empty(self.size)

    @property
    def radial_abscissas(self):
        return np.array(self._radial_abscissas[:self.n_radial])

    @property
    def radial_weights(self):
        return np.array(self._radial_weights[:self.n_radial])

    @property
    def angular_theta(self):
        return np.array(self._angular_theta[:self.n_angular])

    @property
    def angular_phi(self):
        return np.array(self._angular_phi[:self.n_angular])

    @property
    def angular_weights(self):
        return np.array(self._angular_weights[:self.n_angular])

    def free(self):
        """
        Free memory allocated by C routines. Always use it after using a BeckeGrid object.
        """
        napmo_library.grid_free(byref(self))

    def weight_c(self, r, particleID, csystem):
        """Computes the Becke weights :math:`w(r)` at point ``r`` for particle ``particleID`` as described in eq. 22 Becke, 1988
        using C routine.

        References:
            Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

        Args:
            r (array[3]): Point of the grid in which the weight will be calculated.
            particleID (int): The particle index who owns the ``r`` point.
            csys (CBinding) : Information of the system in C structure. See CBinding.

        Returns:
            P (float64): The value of cell_function (eq. 13, Becke, 1988) at point ``r``
        """
        aux = (c_double * 3)()

        for i in range(3):
            aux[i] = c_double(r[i])

        return napmo_library.grid_weights(byref(csystem), aux, particleID)

    def integrate_c(self, csys):
        """
        Perform an integration of function :math:`F(r)` using BeckeGrid. (Function coded on C)

        Args:
            csys (CBinding) : Information of the system in C structure. See CBinding.

        Info:
            So far only implements density integration as a test, more functionals soon.
        """
        return napmo_library.grid_integrate(byref(csys), byref(self))

    def integrate(self, molecule, F, args=()):
        """
        Perform an integration of function :math:`F(r)` using BeckeGrid.

        Args:
            molecule (MolecularSystem): Molecular system to be used in the calculation.
            F (function): Functional to be integrated.
            args (sequence, optional): Extra arguments to pass to F.
        """
        r = np.zeros([3], dtype=np.float64)
        integral = 0.0
        integral_new = 0.0
        system = CBinding(molecule.get('atoms'))

        for i in range(len(molecule.get('atoms'))):
            particle = molecule.get('atoms')[i]
            atgrid = AtomicGrid(
                self.n_radial, self.n_angular, particle.get("origin"), particle.get("symbol"))

            for j in range(atgrid.size):
                r = atgrid.points[j, :]
                atgrid.becke_weights[j] = self.weight_c(
                    r, i, system)

            f = F(atgrid.points, *args)
            integral += atgrid.integrate(f)

        return integral

    def expand(self):
        counter = 0
        for r in range(self.n_radial):
            for a in range(self.n_angular):
                self.xyz[counter, 0] = self._radial_abscissas[r]
                self.xyz[counter, 1] = self._angular_theta[a]
                self.xyz[counter, 2] = self._angular_phi[a]
                self.w[counter] = self._radial_weights[
                    r] * self._angular_weights[a]
                counter += 1
        self.expanded = True

    def convert(self):
        if self.spherical:

            if not self.expanded:
                self.expand()

            aux3 = np.sin(self.xyz[:, 1])
            x = self.xyz[:, 0] * aux3 * np.cos(self.xyz[:, 2])
            y = self.xyz[:, 0] * aux3 * np.sin(self.xyz[:, 2])
            z = self.xyz[:, 0] * np.cos(self.xyz[:, 1])

            self.xyz[:, 0] = x
            self.xyz[:, 1] = y
            self.xyz[:, 2] = z

            self.spherical = False

        else:

            xy = self.xyz[:, 0]**2 + self.xyz[:, 1]**2
            r = np.sqrt(xy + self.xyz[:, 2]**2)
            theta = np.arctan2(np.sqrt(xy), self.xyz[:, 2])
            phi = np.arctan2(self.xyz[:, 1], self.xyz[:, 0])

            self.xyz[:, 0] = r
            self.xyz[:, 1] = theta
            self.xyz[:, 2] = phi

            self.spherical = True

    def move(self, coord=[0.0, 0.0, 0.0], scaling_factor=1.0):

        self.spherical = True
        self.expanded = False
        self.convert()

        if scaling_factor != 1.0:
            self.xyz[:, :] *= scaling_factor

        self.xyz += coord

    def show(self):
        """
        Prints information of the object.
        """
        print("Grid Information:")
        print("-----------------")
        print("Radial Points: ", self.n_radial)
        print("Angular Points: ", self.n_angular)

#################################################
# Interface to napmo_library
#################################################
napmo_library.grid_weights.restype = c_double
napmo_library.grid_weights.argtypes = [
    POINTER(CBinding),
    POINTER(c_double * 3),
    c_int
]
napmo_library.grid_integrate.restype = c_double
napmo_library.grid_integrate.argtypes = [
    POINTER(CBinding),
    POINTER(BeckeGrid)
]
#################################################
