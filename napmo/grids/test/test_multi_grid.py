# file: test_multi_grid.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np


def test_napmo_multi_grid():
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
    system = napmo.NAPMO(data, pprint=False)

    # Create Manager
    mgrid = napmo.MultiGrid(system.system.size_species)

    for p, key in enumerate(system.system.keys()):
        particle = system.system.get(key, {})
        if particle.get('is_electron'):
            key = 'e-'

        aux = system.data.scf.get('grid').get(key, None)

        # Add grids to it
        mgrid.add_grid(napmo.BeckeGrid(particle,
                                       aux.get('nrad', 100),
                                       aux.get('nang', 110),
                                       rtransform=aux.get('rtransform', None)))

    mgrid.show()

    # Sanity checks
    assert mgrid.ngrids == 2
    assert mgrid.get_grid('x') is None
    assert mgrid.get_grid('e-')._symbol == 'e-'
    assert mgrid.get_grid(gid=0)._symbol == 'e-'
    assert mgrid.get_grid_id('e-') == 0
    assert mgrid.get_grid_id('H_1') == 1
    assert mgrid.get_grid_id('x') is None

    # common points
    assert mgrid.get_common_points('e-', 'X') is None
    cpoints = mgrid.get_common_points('e-', 'H_1')
    assert cpoints.shape[0] == 5500
    assert np.allclose(cpoints[:, 0], np.arange(5500, 11000, dtype=np.int))


if __name__ == '__main__':
    test_napmo_multi_grid()
