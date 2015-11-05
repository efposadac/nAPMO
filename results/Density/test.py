
from napmo.interfaces.molecular_system import MolecularSystem
from napmo.grids.becke import GridBecke
from napmo.grids.becke_grid import BeckeGrid

import numpy as np
import os

file_dens = os.path.join(os.path.dirname(__file__), 'H_dens.dat')
P = np.array(np.loadtxt(file_dens), order='F', dtype=np.float64)
os.system('cp ' + file_dens + ' data.dens')

basis_file = os.path.join(os.path.dirname(__file__), "STO-3G.json")

molecule = MolecularSystem()
molecule.add_atom('H', [0.0, 0.0, 0.371],
                  basis_kind='GTO', basis_file=basis_file)
molecule.add_atom('H', [0.0, 0.0, -0.371],
                  basis_kind='GTO', basis_file=basis_file)

molecule.show()

system = CBinding(molecule.get('atoms'))


def rho(coord, molecule, P=P):
    basis = molecule.get_basis_set('e-')
    occupation = molecule.n_occupation('e-')
    bvalue = np.array(basis.compute(coord))
    output = np.empty(coord.shape[0])
    for i in range(coord.shape[0]):
        output[i] = bvalue[:, i].dot(P.dot(bvalue[:, i]))
    return output

angularPoints = 194
radialPoints = 100

grid = []
for i in range(len(molecule.get('atoms'))):
    rm = molecule.get('atoms')[i].get('atomic_radii_2')
    grid.append(BeckeGrid(radialPoints, angularPoints))
    grid[-1].move(molecule.get('atoms')[i].get('origin'), rm)

newgrid = GridBecke(molecule, radialPoints, angularPoints)
newp = newgrid.becke_weights()
p = np.empty(newp.shape)

for i in range(newgrid.ncenter):
    print(np.allclose(grid[i].xyz, newgrid.atgrids[i].points))
    print(np.allclose(grid[i].w, newgrid.atgrids[i].weights))
    for j in range(grid[i].size):
        r = grid[i].xyz[j, :]
        p[j, i] = grid[i].weight_c(r, i, system)
    print(np.allclose(newp[:, i], p[:, i]))
