# file: test_xc_interface.py
# nAPMO package
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

import napmo
import numpy as np
import pathlib


def test_xc_interface_polarized():
    pwd = pathlib.Path(__file__).parent.absolute()
    print(pwd)

    f = napmo.Functional('Li')
    rho = np.fromfile(str(pwd) + '/rho.txt', sep=" ")

    c_zk_ref = np.fromfile(str(pwd) + '/c_zk.txt', sep=" ").reshape(1, 22000)
    c_vrho_ref = np.fromfile(str(pwd) + '/c_vrho.txt', sep=" ").reshape(1, 22000)

    x_zk_ref = np.fromfile(str(pwd) + '/x_zk.txt', sep=" ")
    x_vrho_ref = np.fromfile(str(pwd) + '/x_vrho.txt', sep=" ")

    c_zk, c_vrho = f.compute_correlation(rho)
    assert(np.allclose(c_zk, c_zk_ref))
    assert(np.allclose(c_vrho, c_vrho_ref))

    x_zk, x_vrho = f.compute_exchange(rho)
    assert(np.allclose(x_zk, x_zk_ref))
    assert(np.allclose(x_vrho, x_vrho_ref))


# if __name__ == '__main__':
#     test_xc_interface_polarized()
