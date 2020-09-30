# file: nnuclear.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

from __future__ import division
from __future__ import print_function

import numpy as np
import napmo


def compute_nuclear(grid, point_charges):
    """
    Computes the Nuclear attraction repulsion operator on the grid

    :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|}`

    Args:
        grid (BeckeGrid) : The Molecular grid
        point_charges (list) : List of dict with the information of the point charges.
                               It should have the keys ``atomic_number`` or ``charge`` and ``origin``.

    Return:
        V (ndarray) : :math:`\sum_{A} \dfrac{Z_{A}}{|r_i -R_A|}` for all points in the grid.
    """
    V = np.array([np.array([
        point.get('atomic_number', point.get('charge')) /
        np.sqrt(((point.get('origin') - r)**2).sum())
        for point in point_charges]).sum()
        for r in grid.points])

    return V
