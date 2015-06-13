# file: becke_grid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np

from interfaces.stack import Stack
from utilities.constants import *
from utilities.numerical_integration import *


class BeckeGrid(object):
    """
    This class creates the Becke grid.

    References:
        Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

    Args:
            n_radial (int, optional): Number of radial points. Default is 15
            n_angular (int, optional): Number of angular points. Default is 110
    """

    def __init__(self, n_radial=40, n_angular=110):
        super(BeckeGrid, self).__init__()
        self.n_radial = n_radial
        self.n_angular = n_angular
        self.radial_abscissas, self.radial_weights = chebgauss(n_radial)
        self.angular_theta, self.angular_phi, self.angular_weights = lebedev_q(n_angular)

    def weight(self, r, particleID, particle_stack):
        """Computes the Becke weights :math:`w(r)` at point ``r`` for particle ``particleID`` as described in eq. 22 Becke, 1988
        using Python routine.


        References:
            Becke, A. D. A multicenter numerical integration scheme for polyatomic molecules. J. Chem. Phys. 88, 2547 (1988).

        Args:
            r (array[3]): Point of the grid in which the weight will be calculated.
            particleID (int): The particle index who owns the ``r`` point.
            particle_stack (Stack): stack of particles of same species as ``particleID``

        Returns:
            P (float64): The value of cell_function (eq. 13, Becke, 1988) at point ``r``
        """
        assert isinstance(particleID, int)
        assert isinstance(particle_stack, Stack)

        P = np.ones([len(particle_stack)], dtype=np.float64)

        for i in range(len(particle_stack)):
            P[i] = 1.0
            for j in range(len(particle_stack)):
                if i != j:
                    # Internuclear distance (R_ij eq. 11)
                    x_i = particle_stack[i].get('origin')
                    x_j = particle_stack[j].get('origin')
                    r_i = np.sqrt(np.dot(r - x_i, r - x_i))
                    r_j = np.sqrt(np.dot(r - x_j, r - x_j))
                    R_ij = np.sqrt(np.dot(x_i - x_j, x_i - x_j))

                    # \mu_ij eq. 11
                    mu_ij = (r_i - r_j) / R_ij

                    # Atomic size adustment. see apendix, Becke, 1988.
                    rm_i = particle_stack[i].get('atomic_radii') * ANGSTROM_TO_BOHR
                    rm_j = particle_stack[j].get('atomic_radii') * ANGSTROM_TO_BOHR

                    # eq. A4
                    chi = rm_i / rm_j
                    # eq. A6
                    u_ij = (chi - 1.0) / (chi + 1.0)
                    # eq. A5
                    a_ij = u_ij / ((u_ij * u_ij) - 1.0)
                    # eq. A3
                    if (np.abs(a_ij) > 0.50):
                        a_ij = 0.50 * a_ij / np.abs(a_ij)
                    # eq. A2
                    nu_ij = mu_ij + a_ij * (1.0 - mu_ij * mu_ij)

                    P[i] = P[i] * self.step_function(3, nu_ij)

        # eq. 22
        return P[particleID] / np.sum(P)

    def integrate(self, particle_stack, F):
        """
        Perform an integration of function :math:`F(r)` using BeckeGrid. (Function coded on Python)
        """
        r = np.zeros([3], dtype=np.float64)
        integral = 0.0
        for n in range(self.n_radial):
            x = self.radial_abscissas[n]
            r_w = self.radial_weights[n]
            for m in range(int(self.n_angular)):
                phi = self.angular_phi[m]
                theta = self.angular_theta[m]
                a_w = self.angular_weights[m]
                for i in range(len(particle_stack)):
                    particle = particle_stack[i]
                    rm = particle.get('atomic_radii') * ANGSTROM_TO_BOHR

                    if particle.get("atomic_number") != 1:
                        rm *= 0.5

                    aux = (x + 1.0) * 0.5
                    aux2 = aux * aux
                    rad = -rm * np.log(1.0 - (aux2 * aux2))

                    aux3 = np.sin(theta)
                    r[0] = rad * aux3 * np.cos(phi)
                    r[1] = rad * aux3 * np.sin(phi)
                    r[2] = rad * np.cos(theta)

                    # Move grid points
                    r += particle.get("origin")

                    # Calculate Becke weigths
                    p = self.weight(r, i, particle_stack)

                    # This factor comes from the variable change r ! --> x ,
                    # and from using chebyshev-gauss radial quadrature of second order
                    factor = 2.0 * rm * (aux2 * aux)/(np.sqrt(1.0 - x * x) * (1.0 - (aux2 * aux2)))
                    integral += 4.0 * np.pi * rad * rad * p * r_w * a_w * F(r, particle_stack) * factor

        return integral

    def step_function(self, order, mu):
        """
        Iterated cutoff profile. eq. 21, Becke 1988.
        """
        # eq. 19
        def P(_mu):
            return 1.5 * _mu - 0.5 * _mu * _mu * _mu

        f = mu
        for k in range(order):
            f = P(f)
        return 0.5 * (1. - f)

    def show(self):
        """
        Prints information of the object.
        """
        print("Grid Information:")
        print("-----------------")
        print("Radial Points: ", self.get_n_radial())
        print("Angular Points: ", self.get_n_angular())
