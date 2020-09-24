# file: test_multigrid.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# fernando.posada@temple.edu

import napmo
import numpy as np


def test_napmo_multigrid():
    file = ("""
# Molecule definition:

molecule [e-:0:1] {
    H    0.000000  0.00000  0.367
    H_1  0.000000  0.00000 -0.367
}

basis {
   e- STO-3G
   H_1 DZSNB
}

scf {
    maxiter = 100
    dft
}

grid {
    e- [50, 110]
    H_1 [50, 110] 
}
""")

    data = napmo.InputParser(file)
    system = napmo.NAPMO(data)
    solver = napmo.DFT(system.system, options=system.data.scf)

    # Create manager
    mgrid = napmo.MultiGrid(system.system.size_species)

    # Add grids to it
    for psi in solver.PSI:
        mgrid.add_grid(psi._grid)

    mgrid.show()

    # Sanity checks
    assert mgrid.ngrids == 2
    assert mgrid.get_grid('x') == None
    assert mgrid.get_grid('e-')._symbol == 'e-'
    assert mgrid.get_grid_id('e-') == 0
    assert mgrid.get_grid_id('H_1') == 1
    assert mgrid.get_grid_id('x') == None

    # common points
    assert mgrid.get_common_points('e-', 'X') == None
    cpoints = mgrid.get_common_points('e-', 'H_1')
    assert cpoints.shape[0] == 5500
    assert np.allclose(cpoints, np.arange(5500, 11000, dtype=np.int))


# if __name__ == '__main__':
#     test_napmo_multigrid()