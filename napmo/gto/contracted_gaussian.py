# file: contracted_gaussian.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
from ctypes import *
import scipy.misc

import napmo


class ContractedGaussian(object):

    """
    Defines a linear combination of Cartesian Gaussian type orbitals GTO.

    A contracted Gaussian function is just a linear combination of primitive
    Gaussians (also termed primitives) centered at the same center
    :math:`{\\bf A}` and with the same momentum indices  :math:`{\\bf n}`
    but with different exponents :math:`\zeta_i`:

    :math:`\phi ({\\bf r}; {\\bf \zeta}, {\\bf C}, {\\bf n}, {\\bf A}) =
    (x - A_x)^{n_x} (y - A_y)^{n_y} (z - A_z)^{n_z} \\times \sum_{i=1}^M C_i
    \exp[-\zeta_i ({\\bf r}-{\\bf A})^2]`

    Contracted Gaussians form shells the same way as primitives. The
    contraction coefficients :math:`\\bf C` already include normalization
    constants so that the resulting combination is properly normalized.
    Published contraction coefficients :math:`\\bf c` are linear coefficients
    for normalized primitives, hence the normalization-including contraction
    coefficients :math:`\\bf C` have to be computed from them as;

    :math:`C_i = c_i N(\zeta_i,{\\bf n})`

    where :math:`N` is:

    :math:`N = \\dfrac{1}{\sqrt{<\phi | \phi>}}`

    Args:
        exponents (ndarray): GTO exponent.
        coefficients (ndarray): GTO coefficients.
        origin (ndarray) : coordinates (cartesian)
        l (ndarray) : :math:`\\bf n`. Angular moment (x, y, and z components)
    """

    def __init__(self,
                 exponents=np.array([0.5], dtype=np.float64),
                 coefficients=np.array([1.0], dtype=np.float64),
                 origin=np.zeros(3, dtype=np.float64),
                 l=np.zeros(3, dtype=np.int32)):

        super(ContractedGaussian, self).__init__()

        assert isinstance(exponents, np.ndarray) or isinstance(exponents, list)

        self._primitives = [
            napmo.PrimitiveGaussian(exponent, coefficient, l, origin)
            for (exponent, coefficient) in zip(exponents, coefficients)]

        aux = np.array([p._this for p in self._primitives], dtype=c_void_p)

        self._this = napmo.cext.ContractedGaussian_new(aux, len(exponents))

    def overlap(self, other):
        """
        Calculates the overlap integral between two contractions.

        Args:
            other (ContractedGaussian) : Contracted function to perform :math:`<\phi_{self} | \phi_{other}>`

        """

        return napmo.cext.ContractedGaussian_overlap(self._this, other._this)

    def compute(self, coord):
        """
        Computes the value of the contracted Gaussian at ``coord``.

        Args:
            coord (ndarray) : array with the points where the function will be
            calculated.
        """
        n_coord = coord.shape[0]
        output = np.empty(n_coord)
        napmo.cext.ContractedGaussian_compute(
            self._this, coord, output, n_coord)

        return output

    @property
    def l(self):
        """
        The angular moment of the object
        """
        output = np.zeros(3, dtype=np.int32)
        napmo.cext.ContractedGaussian_get_l(self._this, output)
        return output

    @property
    def origin(self):
        """
        The center of the function
        """
        output = np.zeros(3)
        napmo.cext.ContractedGaussian_get_origin(self._this, output)
        return output

    @property
    def nprim(self):
        return napmo.cext.ContractedGaussian_get_nprim(self._this)

    @property
    def normalization(self):
        return napmo.cext.ContractedGaussian_get_norma(self._this)

    def _show_compact(self):
        """
        contractions information of the object in a "compact" way
        """

        out = ''
        out += "".join([p._show_compact()
                        for p in self._primitives if p.l[0] == sum(p.l)])

        return out

    def __repr__(self):

        out = """
==================================================
Object: {0:9s}
--------------------------------------------------
Origin: {1:<18.14f} {2:<18.14f} {3:<18.14f}
l:      {4:<3d} {5:<3d} {6:<3d}
nprim: {7:<10d}
norma:  {8:<18.14f}
--------------------------------------------------
""".format(
            type(self).__name__,
            self.origin[0],
            self.origin[1],
            self.origin[2],
            self.l[0],
            self.l[1],
            self.l[2],
            self.nprim,
            self.normalization)

        out += """
  {0:<3s} {1:>10s} {2:>10s} {3:>10s}
  ------------------------------------
""".format(
            "l",
            "zeta",
            "Coeff",
            "Norma")

        out += "".join([p._show_compact()
                        for p in self._primitives])

        return out
