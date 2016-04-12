import numpy as np
import matplotlib.pyplot as plt
import time

from napmo.utilities.lebedev import lebedev_get_order
from napmo.system.molecular_system import MolecularSystem
from napmo.grids.becke import BeckeGrid
from napmo.grids.poisson_solver import poisson_solver
from napmo.utilities.cubic_spline import CubicSpline as CS

from horton import *

basis_file = "/home/fernando/PhD/Develop/nAPMO/results/Coulomb/TEST.json"
molecule = MolecularSystem()
molecule.add_atom("He", [0.0, 0.0, 0.0],
                  basis_kind="GTO", basis_file=basis_file)
molecule.show()

# Grid definition
angularPoints = 110
radialPoints = 100

# Grids
grid = BeckeGrid(molecule, radialPoints, angularPoints)

rmin = grid.atgrids[-1].radial_grid.points.min()
rmax = grid.atgrids[-1].radial_grid.points.max()

print(rmin, rmax)

h_grid_spec = AtomicGridSpec(
    'power:' + str(rmin) + ':' + str(rmax) + ':' + str(radialPoints) + ':' + str(angularPoints))

molgrid = BeckeMolGrid(np.array([[0.0, 0.0, 0.0]]), np.array([2]), np.array(
    [2.0]), agspec=h_grid_spec, random_rotate=False, mode='keep')

hatgrid = molgrid.subgrids[-1]
natgrid = grid.atgrids[-1]

# Basis and functional
basis = molecule.get('atoms')[-1].get('basis')
a = basis.get('function')[0]
b = basis.get('function')[0]


def p_ab(coord):
    return a.compute(coord) * b.compute(coord)

hrho = p_ab(molgrid.points)
nrho = p_ab(grid.points)

# tests
h = hatgrid.get_spherical_average(hrho)
n = natgrid.spherical_average(nrho)

np.allclose(h, n)
