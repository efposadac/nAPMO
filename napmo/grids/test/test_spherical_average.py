# file: test_spherical_average.py
# nAPMO package
# Copyright (c) 2016, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu


import numpy as np
import napmo
import matplotlib.pyplot as plt


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
        print(sa_fn.shape)
        print(abs(sa_fn - sa_check).max())
        assert abs(sa_fn - sa_check).max() < 1e-10


def test_spherical_decomposition_hydrogen_1s():
    rtf = napmo.PowerRadialTransform(1e-3, 2e1, 100)
    origin = np.random.uniform(-1, 1, 3)
    grid = napmo.AtomicGrid(100, 110, origin, 'H', rtransform=rtf)

    distances = np.sqrt(((origin - grid.points)**2).sum(axis=1))
    fn = np.exp(-2 * distances) / np.pi

    sa_fns = grid.spherical_expansion(4, fn)
    sa_check = np.exp(-2 * grid.radial_grid.points) / \
        np.pi * np.sqrt(4 * np.pi)

    # plt.plot(grid.radial_grid.points, sa_fns[:, 0], label='napmo')
    # plt.plot(grid.radial_grid.points, sa_check, label='patron')
    # plt.legend()
    # plt.show()

    assert abs(sa_fns[:, 0] - sa_check).max() < 1e-10

    for sa_fn in sa_fns[:, 1]:
        assert abs(sa_fn).max() < 1e-10


def test_spherical_decomposition_hydrogen_1pz():
    # density of the 1pz orbital
    origin = np.random.uniform(-1, 1, 3)
    rtf = napmo.PowerRadialTransform(1e-3, 2e1, 100)
    grid = napmo.AtomicGrid(100, 110, origin, 'H', rtransform=rtf)

    z = grid.points[:, 2] - origin[2]
    distances = np.sqrt(((origin - grid.points)**2).sum(axis=1))
    fn = np.exp(-distances) / (32.0 * np.pi) * z**2

    sa_fns = grid.spherical_expansion(4, fn)

    # s
    sa_check = np.exp(-grid.radial_grid.points) / (32.0 * np.pi) * \
        (1.0 / 3.0) * grid.radial_grid.points**2 * np.sqrt(4 * np.pi)
    assert abs(sa_fns[:, 0] - sa_check).max() < 1e-10
    # p
    for sa_fn in sa_fns[4:1]:
        assert abs(sa_fn).max() < 1e-10
    # d
    sa_check = np.exp(-grid.radial_grid.points) / (32.0 * np.pi) * \
        (2.0 / 15.0) * grid.radial_grid.points**2 * np.sqrt(5 * 4 * np.pi)

    assert abs(sa_fns[:, 4] - sa_check).max() < 1e-10  # 6 or 4???

# test_spherical_average_H1s()
# test_spherical_decomposition_hydrogen_1s()
# test_spherical_decomposition_hydrogen_1pz()
