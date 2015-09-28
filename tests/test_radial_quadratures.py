# file: test_radial_quadratures.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

import napmo.utilities.radial_quadratures as nint
import numpy as np


def test_radial_quadrature_chebgauss():
    _r = np.array([8.66025404e-01, 5.00000000e-01, 6.12323400e-17, -5.00000000e-01, -8.66025404e-01])
    _w = np.array([0.13089969, 0.39269908, 0.52359878, 0.39269908, 0.13089969])
    r, w = nint.chebgauss(5)
    np.testing.assert_allclose(_r, r)
    np.testing.assert_allclose(_w, w)


def test_radial_quadrature_chebgauss_transformed():

    try:
        nint.chebgauss_transformed(-1)
        assert False, 'Expecting Failure!'
    except:
        assert True

    # test for 15 points.
    _rq = np.zeros([15, 2])

    _rq[:, 0] = np.array([
        0.99990271322752600, 0.99705379102173697, 0.97959347194838908, 0.92441318157838759,
        0.80461983162814810, 0.60315708641633381, 0.32492907290728490, 0.00000000000000000,
        -0.32492907290728490, -0.60315708641633381, -0.80461983162814810, -0.92441318157838759,
        -0.97959347194838908, -0.99705379102173697, -0.99990271322752600], dtype=np.float64)

    _rq[:, 1] = np.array([
        4.8286046422502072e-004, 7.1488698022420796e-003, 3.1756645389712214e-002, 8.3333333333333301e-002,
        0.15931778951140874, 0.24285113019775792, 0.30844270463465390, 0.33333333333333331,
        0.30844270463465390, 0.24285113019775792, 0.15931778951140874, 8.3333333333333301e-002,
        3.1756645389712214e-002, 7.1488698022420796e-003, 4.8286046422502072e-004], dtype=np.float64)

    rq = nint.chebgauss_transformed(15)

    np.testing.assert_allclose(rq, _rq)


def test_radial_quadrature_chebgauss_transformed_integrate():
    eps = 1.0e-12

    f = lambda x: np.exp(-6.793 * x * x) / (1.000001 - (x * x))
    i, err = nint.chebgauss_transformed_integrate(f, eps)
    np.testing.assert_allclose(i, 0.762416793289, rtol=eps)
    assert err < eps

    f = lambda x: 1 / (x**4 + x**2 + 0.9)
    i, err = nint.chebgauss_transformed_integrate(f, eps)
    np.testing.assert_allclose(i, 1.582232963730, rtol=eps)
    assert err < eps

    try:
        i, err = nint.chebgauss_transformed_integrate(f, 1.0e-20)
        assert False, 'Expecting Failure!'
    except:
        assert True

    try:
        i, err = nint.chebgauss_transformed_integrate(f, max_iter=5)
        assert False, 'Expecting Failure!'
    except:
        assert True


# test_radial_quadrature_chebgauss_transformed()
# test_radial_quadrature_chebgauss_transformed_integrate()
# test_radial_quadrature_chebgauss()
