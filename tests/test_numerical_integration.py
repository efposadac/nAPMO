# file: test_numerical_integration.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

from utilities import numerical_integration as nint
import numpy as np


def test_numerical_integration_chebgauss_t_q_w():

    try:
        nint.chebgauss_t_q_w(-1)
        assert False, 'Assert exception expected'
    except:
        pass

    # test for 15 points.
    _q = [
        0.99990271322752600, 0.99705379102173697, 0.97959347194838908, 0.92441318157838759,
        0.80461983162814810, 0.60315708641633381, 0.32492907290728490, 6.4969530541989581e-017,
        -0.32492907290728490, -0.60315708641633381, -0.80461983162814810, -0.92441318157838759,
        -0.97959347194838908, -0.99705379102173697, -0.99990271322752600]

    _w = [
        4.8286046422502072e-004, 7.1488698022420796e-003, 3.1756645389712214e-002, 8.3333333333333301e-002,
        0.15931778951140874, 0.24285113019775792, 0.30844270463465390, 0.33333333333333331,
        0.30844270463465390, 0.24285113019775792, 0.15931778951140874, 8.3333333333333301e-002,
        3.1756645389712214e-002, 7.1488698022420796e-003, 4.8286046422502072e-004]

    q, w = nint.chebgauss_t_q_w(len(_w))

    for i in xrange(len(_w)):
        assert np.abs(_q[i] - q[i]) < 1.0e-12
        assert np.abs(_w[i] - w[i]) < 1.0e-12

    # Data from FORTRAN code.  The difference 10-7 is due a the difference of Ln(2.0)
    # Between Fortran and the rest of the world. :-(

    _q = [
        14.3273967744679, 9.40692449540646, 6.61482542993771, 4.72572151461268,
        3.35564405121858, 2.33336004768934, 1.56688900481176, 1.00000000000000,
        0.594084870800892, 0.319084206016552, 0.148305057079304, 5.558141605618139E-002,
        1.479581291802201E-002, 2.126810172847216E-003, 7.018202663133799E-005]

    _w = [
        7.16048418042371, 3.50064746235206, 2.24512247145824, 1.59054963445867,
        1.17640897813369, 0.882868531238085, 0.659173342078948, 0.480898345641558,
        0.335858551413324, 0.218543849019213, 0.127365875179677, 6.247337503258128E-002,
        2.314371882021904E-002, 5.164427232597799E-003, 3.483271414439691E-004]

    q, w = nint.chebgauss_t_q_w(len(_w), rescale=True)

    for i in xrange(len(_w)):
        print _q[i], q[i], np.abs(_q[i] - q[i])
        assert np.abs(_q[i] - q[i]) < 1.0e-7
        assert np.abs(_w[i] - w[i]) < 1.0e-7


def test_numerical_integration_chebgauss_integrate():
    eps = 1.0e-12

    f = lambda x: np.exp(-6.793 * x * x) / (1.000001 - (x*x))
    i, err = nint.chebgauss_integrate(f, eps)
    assert np.abs(i-0.762416793289) < eps
    assert err < eps

    f = lambda x: 1 / (x**4 + x**2 + 0.9)
    i, err = nint.chebgauss_integrate(f, eps)
    assert np.abs(i-1.582232963730) < eps
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


def test_numerical_integration_lebedev_q_w():
    # test for 14 points.
    _x = [
        1.0000000000000000, -1.0000000000000000, 0.0000000000000000, 0.0000000000000000,
        0.0000000000000000, 0.0000000000000000, 0.57735026918962573, -0.57735026918962573,
        0.57735026918962573, -0.57735026918962573, 0.57735026918962573, -0.57735026918962573,
        0.57735026918962573, -0.57735026918962573]

    _y = [
        0.0000000000000000, 0.0000000000000000, 1.0000000000000000, -1.0000000000000000,
        0.0000000000000000, 0.0000000000000000, 0.57735026918962573, 0.57735026918962573,
        -0.57735026918962573, -0.57735026918962573, 0.57735026918962573, 0.57735026918962573,
        -0.57735026918962573, -0.57735026918962573]

    _z = [
        0.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000,
        1.0000000000000000, -1.0000000000000000, 0.57735026918962573, 0.57735026918962573,
        0.57735026918962573, 0.57735026918962573, -0.57735026918962573, -0.57735026918962573,
        -0.57735026918962573, -0.57735026918962573]

    _w = [
        6.6666666666666666E-002, 6.6666666666666666E-002, 6.6666666666666666E-002, 6.6666666666666666E-002,
        6.6666666666666666E-002, 6.6666666666666666E-002, 7.4999999999999997E-002, 7.4999999999999997E-002,
        7.4999999999999997E-002, 7.4999999999999997E-002, 7.4999999999999997E-002, 7.4999999999999997E-002,
        7.4999999999999997E-002, 7.4999999999999997E-002]

    x, y, z, w = nint.lebedev_q_w(len(_x))

    for i in xrange(len(_x)):
        assert np.abs(_x[i] - x[i]) < 1.0e-12
        assert np.abs(_y[i] - y[i]) < 1.0e-12
        assert np.abs(_z[i] - z[i]) < 1.0e-12
        assert np.abs(_w[i] - w[i]) < 1.0e-12

    # Number of points supported:
    # 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770
    try:
        x, y, z, w = nint.lebedev_q_w(15)
        assert False, 'Number of points not supported!'
    except:
        pass
