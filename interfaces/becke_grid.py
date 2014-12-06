# file: becke_grid.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from __future__ import division
import numpy as np
from copy import deepcopy

from utilities import lebedev

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
    	""" Computes the radial distribution for a Becke-like grid 
    	"""
        r, w = self.chebgauss(n)

        # domain transformation ... (check)
        w = (w/(1.-r))/np.log(2.)
    	r = 1.- np.log(1.-r)/np.log(2.0)

    	return r, w

    def chebgauss(self, n):
		"""Computes the abscissas and weights for transformed Gauss-Chebyshev quadrature of second kind.
		
		Special case for the integral I f(x) dx. since Gauss-Chebyshev is for I f(x) sqrt(1-x^2). 

		see: Perez-Jorda, J., San-Fabian, E. & Moscardo, F. A simple, reliable and efficient scheme for 
		automatic numerical integration. Comput. Phys. Commun. 70, 271-284 (1992).

		returns the weights and abscissas to perform integration
		"""
		w = np.zeros(n, dtype=np.float64)
		r = np.zeros(n, dtype=np.float64)

		count = 1
		for i in xrange(int((n+1)/2)):
			d = n + 1.
			s = np.sin((count*np.pi)/d)
			c = np.cos((count*np.pi)/d)
			s2 = s*s
			s4 = s2 * s2
			r[i] = 1. + (2./np.pi) * (1. + (2./3.) * s2) * c * s - (2. * count)/d
			r[n-i-1] = -r[i]
			w[i] = 16./(3.*d) * s4
			w[n-i-1] = w[i]
			count += 1

		return r, w

    def angular_distr(self, n):
		""" Computes the angular distribution for a Becke-like grid 
		param x, y, z: Cartesian coordinates
		param w: weights
		see: V.I. Lebedev, and D.N. Laikov, Doklady Mathematics, 59, No. 3, 477 (1999)
		"""

		x = np.zeros(n, dtype=np.float64)
		y = np.zeros(n, dtype=np.float64)
		z = np.zeros(n, dtype=np.float64)
		w = np.zeros(n, dtype=np.float64)

		if lebedev.lebedev_compute(x, y, z, w, n):
			return x, y, z, w
		else:
			print "ERROR!!! Invalid number of angular points."
			print "Info: Only supported:..."
			print "Info: 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770"
			return 0, 0, 0, 0

