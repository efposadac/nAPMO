# file: test_spherical_average.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co


import numpy as np
import napmo


def test_spherical_average_H1s():
    # density of the 1s orbital

    rtf = napmo.PowerRadialTransform(1e-3, 2e1, 100)
    origin = np.random.uniform(-1, 1, 3)
    grid = napmo.AtomicGrid(100, 110, origin, 'H', rtransform=rtf)

    distances = np.sqrt(((origin - grid.points)**2).sum(axis=1))
    fn = np.exp(-2 * distances) / np.pi

    x = grid.points[:, 0] - grid.origin[0]
    y = grid.points[:, 1] - grid.origin[1]
    z = grid.points[:, 2] - grid.origin[2]

    sa_check = np.exp(-2 * grid.radial_grid.points) / np.pi

    for cx, cy, cz, cxxx in (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (1, 1, 0, 0), (0, 1, 0, 1):
        sa_fn = grid.spherical_average(
            fn + cx * x + cy * y + cz * z + cxxx * x * x * x)
        print(abs(sa_fn - sa_check).max())
        assert abs(sa_fn - sa_check).max() < 1e-10


# test_spherical_average_H1s()
