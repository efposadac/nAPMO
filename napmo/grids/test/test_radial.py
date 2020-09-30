# file: test_radial.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu.co

from napmo.grids.radial import RadialGrid
from napmo.grids.radial_transform import *
from napmo.utilities.int1d import StubIntegrator1D


def test_integrate_gauss():
    rtf = PowerRadialTransform(0.0005, 1e1, 100)
    grid = RadialGrid(rtransform=rtf)
    y = np.exp(-0.5 * grid.points**2)

    print(grid.integrate(y, 1))
    print((2 * np.pi)**1.5)

    assert abs(grid.integrate(y, 1) - (2 * np.pi)**1.5) < 1e-9

    grid = RadialGrid(100)
    y = np.exp(-0.5 * grid.points**2)

    print(grid.integrate(y, 1))
    print((2 * np.pi)**1.5)

    assert abs(grid.integrate(y, 1) - (2 * np.pi)**1.5) < 1e-9

# test_integrate_gauss()
