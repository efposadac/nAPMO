# file: lebedev.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

from __future__ import division
from __future__ import print_function

import numpy as np


def lebedev_get_order(n):
    """
    Return the order of the Lebedev quadrature from the number of points ``n``.

    Args:
        n (int) : number of points.

    Returns:
        int: order of lebedev quadrature.
    """
    lebpoints = {
        6: 3, 14: 5, 50: 11, 110: 17, 146: 19, 170: 21, 194: 23, 230: 25, 266: 27, 302: 29, 350: 31, 434: 35,
        590: 41, 770: 47, 974: 53, 1202: 59, 1454: 65, 1730: 71, 2030: 77, 2354: 83,
        2702: 89, 3074: 95, 3470: 101, 3890: 107, 4334: 113, 4802: 119, 5294: 125, 5810: 131
    }

    return lebpoints[n]
