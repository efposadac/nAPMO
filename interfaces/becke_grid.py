# file: becke_grid.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it 

from __future__ import division
import numpy as np
from copy import deepcopy

from utilities import numerical_integration as nint

class BeckeGrid(object):
    """This class creates the Becke grid.
    see: Becke, A. D. A multicenter numerical integration scheme for polyatomic 
    molecules. J. Chem. Phys. 88, 2547 (1988).
    """
    def __init__(self, n_radial=15, n_angular=110):
        super(BeckeGrid, self).__init__()

        self.radial = n_radial
        self.angular = n_angular
        self.size = self.radial * self.angular

        self.x = np.zeros(self.size)
        self.y = np.zeros(self.size)
        self.z = np.zeros(self.size)
        self.w = np.zeros(self.size)

    def radial_distr(self, n):
    	""" Computes the radial distribution for a Becke-like grid (uses Chebyshev-Gauss method.)

    	n : number of radial points.
    	"""
        r, w = nint.chebgauss_transformed(n)

        # domain transformation ... (check)
        w = (w/(1.-r))/np.log(2.)
    	r = 1.- np.log(1.-r)/np.log(2.0)

    	return r, w

    def angular_distr(self, n):
		""" Computes the angular distribution for a Becke-like grid 
		see: V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59, No. 3, 477 (1999)
		x, y, z: Cartesian coordinates
		w: weights
		n : number of angular points
		"""

		x, y, z, w = nint.lebedev_weights(n)
		return x, y, z, w
