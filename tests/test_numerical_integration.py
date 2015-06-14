# file: test_numerical_integration.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it
import os
import sys

lib_path = os.path.abspath('../')
sys.path.append(lib_path)

from utilities import numerical_integration as nint
import numpy as np


def test_numerical_integration_chebgauss():
    _r = np.array([8.66025404e-01, 5.00000000e-01, 6.12323400e-17, -5.00000000e-01, -8.66025404e-01])
    _w = np.array([0.13089969, 0.39269908, 0.52359878, 0.39269908, 0.13089969])
    r, w = nint.chebgauss(5)
    np.testing.assert_allclose(_r, r)
    np.testing.assert_allclose(_w, w)


def test_numerical_integration_chebgauss_rq():

    try:
        nint.chebgauss_rq(-1)
        assert False, 'Assert exception expected'
    except:
        pass

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

    rq = nint.chebgauss_rq(15)

    np.testing.assert_allclose(rq, _rq)

    _rq[:, 0] = [
        14.3273967744679, 9.40692449540646, 6.61482542993771, 4.72572151461268,
        3.35564405121858, 2.33336004768934, 1.56688900481176, 1.00000000000000,
        0.594084870800892, 0.319084206016552, 0.148305057079304, 5.558141605618139E-002,
        1.479581091802201E-002, 2.126807472847216E-003, 7.01792783133799E-005]

    _rq[:, 1] = [
        7.16048418042371, 3.50064746235206, 2.24512247145824, 1.59054963445867,
        1.17640897813369, 0.882868531238085, 0.659173342078948, 0.480898345641558,
        0.335858551413324, 0.218543849019213, 0.127365875179677, 6.247337503258128E-002,
        2.314371882021904E-002, 5.164427232597799E-003, 3.483271414439691E-004]

    rq = nint.chebgauss_rq(15, rescale=True)

    np.testing.assert_allclose(rq, _rq)


def test_numerical_integration_chebgauss_integrate():
    eps = 1.0e-12

    f = lambda x: np.exp(-6.793 * x * x) / (1.000001 - (x * x))
    i, err = nint.chebgauss_integrate(f, eps)
    np.testing.assert_allclose(i, 0.762416793289, rtol=eps)
    assert err < eps

    f = lambda x: 1 / (x**4 + x**2 + 0.9)
    i, err = nint.chebgauss_integrate(f, eps)
    np.testing.assert_allclose(i, 1.582232963730, rtol=eps)
    assert err < eps

    try:
        i, err = nint.chebgauss_integrate(f, 1.0e-20)
        assert False, 'convergence failed expected'
    except:
        pass

    try:
        i, err = nint.chebgauss_integrate(f, max_iter=5)
        assert False, 'convergence failed expected'
    except:
        pass

test_numerical_integration_chebgauss_rq()
test_numerical_integration_chebgauss_integrate()
test_numerical_integration_chebgauss()
