# file: test_radial.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

from napmo.grids.radial import RadialGrid
from napmo.grids.radial_transform import *


def test_integrate_gauss():
    rtf = PowerRadialTransform(0.0005, 1e1, 100)
    grid = RadialGrid(rtransform=rtf)
    assert isinstance(grid.int1d, StubIntegrator1D)

    y = np.exp(-0.5 * grid.points**2)
    assert abs(grid.integrate(y, 1) - (2 * np.pi)**1.5) < 1e-9

# test_integrate_gauss()
