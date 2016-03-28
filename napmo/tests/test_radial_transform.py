# file: test_radial_transform.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from napmo.grids.radial_transform import *
from napmo.utilities.int1d import *

import numpy as np
from nose.tools import assert_raises


def check_consistency(rtf):
    ts = np.random.uniform(0, rtf.size - 1, 200)
    # consistency between radius and radius_array
    rs = rtf.radius(ts)
    print(ts.shape)
    print(rs.shape)
    for i in range(ts.shape[0]):
        assert rs[i] == rtf.radius(ts[i])
    # consistency between deriv and deriv_array
    ds = rtf.deriv(ts)
    for i in range(ts.shape[0]):
        assert ds[i] == rtf.deriv(ts[i])
    # consistency between deriv2 and deriv2_array
    d2s = rtf.deriv2(ts)
    for i in range(ts.shape[0]):
        assert d2s[i] == rtf.deriv2(ts[i])
    # consistency between deriv3 and deriv3_array
    d3s = rtf.deriv3(ts)
    for i in range(ts.shape[0]):
        assert d3s[i] == rtf.deriv3(ts[i])
    # consistency between inv and inv_array
    ts = rtf.inv(rs)
    for i in range(ts.shape[0]):
        assert ts[i] == rtf.inv(rs[i])
    # consistency between inv and radius
    for i in range(ts.shape[0]):
        assert abs(ts[i] - rtf.inv(rtf.radius(ts[i]))) < 1e-10

    ts = np.arange(rtf.size, dtype=float)
    # consistency of get_radius
    assert (rtf.radius_all() == rtf.radius(ts)).all()
    # consistency of deriv_all
    assert (rtf.deriv_all() == rtf.deriv(ts)).all()
    # consistency of deriv2_all
    assert (rtf.deriv2_all() == rtf.deriv2(ts)).all()
    # consistency of deriv3_all
    assert (rtf.deriv3_all() == rtf.deriv3(ts)).all()

    # radius must increase strictly
    radius = rtf.radius_all()
    assert (radius[1:] > radius[:-1]).all()


def check_deriv(rtf):
    ts = np.random.uniform(0, rtf.size - 1, 200)
    eps = 1e-5
    ts0 = ts - eps / 2
    ts1 = ts + eps / 2
    fns = [(rtf.radius, rtf.deriv),
           (rtf.deriv, rtf.deriv2),
           (rtf.deriv2, rtf.deriv3)]
    for fnr, fnd in fns:
        ds = fnd(ts)
        dns = (fnr(ts1) - fnr(ts0)) / eps
        assert abs(ds - dns).max() < 1e-5


def check_chop(rtf1):
    assert rtf1.size == 100
    rtf2 = rtf1.chop(50)
    assert rtf1.__class__ == rtf2.__class__
    assert rtf2.size == 50
    assert abs(rtf1.get_radius()[:50] - rtf2.get_radius()).max() < 1e-8


def check_half(rtf1):
    radius1 = rtf1.get_radius()
    rtf2 = rtf1.half()
    radius2 = rtf2.get_radius()
    assert abs(radius1[1::2] - radius2).max() < 1e-9


def test_identity_basics():
    rtf = IndentityRadialTransform(100)
    assert rtf.radius(0.0) == 0.0
    assert rtf.radius(99.0) == 99.0
    check_consistency(rtf)
    check_deriv(rtf)
    # check_chop(rtf)


# def test_linear_basics():
#     rtf = LinearRTransform(-0.7, 0.8, 100)
#     assert abs(rtf.radius(0) - -0.7) < 1e-15
#     assert abs(rtf.radius(99) - 0.8) < 1e-10
#     check_consistency(rtf)
#     check_deriv(rtf)
#     check_chop(rtf)
#     check_half(rtf)


# def test_exp_basics():
#     rtf = ExpRTransform(0.1, 1e1, 100)
#     assert abs(rtf.radius(0) - 0.1) < 1e-15
#     assert abs(rtf.radius(99) - 1e1) < 1e-10
#     check_consistency(rtf)
#     check_deriv(rtf)
#     check_chop(rtf)
#     check_half(rtf)


# def test_shifted_exp_basics():
#     rtf = ShiftedExpRTransform(0.1, 0.1, 10.0, 100)
#     assert abs(rtf.radius(0) - 0.1) < 1e-15
#     assert abs(rtf.radius(99) - 10.0) < 1e-10
#     check_consistency(rtf)
#     check_deriv(rtf)
#     check_chop(rtf)


def get_power_cases():
    return [
        (1e-3, 1e2),
        (1e-3, 1e3),
        (1e-3, 1e4),
        (1e-3, 1e5),
    ]


def test_power_basics():
    cases = get_power_cases()
    for rmin, rmax in cases:
        rtf = PowerRadialTransform(rmin, rmax, 100)
        assert abs(rtf.radius(99) - rmax) < 1e-9
        check_consistency(rtf)
        check_deriv(rtf)
        # check_chop(rtf)
        # check_half(rtf)


def test_identity_properties():
    rtf = IndentityRadialTransform(100)
    assert rtf.size == 100
    assert isinstance(rtf.get_default_int1d(), SimpsonIntegrator1D)


# def test_linear_properties():
#     rtf = LinearRTransform(-0.7, 0.8, 100)
#     assert rtf.rmin == -0.7
#     assert rtf.rmax == 0.8
#     assert rtf.size == 100
#     assert rtf.alpha > 0
#     assert isinstance(rtf.get_default_int1d(), SimpsonIntegrator1D)


# def test_exp_properties():
#     rtf = ExpRTransform(0.1, 1e1, 100)
#     assert rtf.rmin == 0.1
#     assert rtf.rmax == 1e1
#     assert rtf.size == 100
#     assert rtf.alpha > 0
#     assert isinstance(rtf.get_default_int1d(), SimpsonIntegrator1D)


# def test_shifted_exp_properties():
#     rtf = ShiftedExpRTransform(0.2, 0.1, 1e1, 100)
#     assert rtf.rmin == 0.2
#     assert rtf.rshift == 0.1
#     assert rtf.rmax == 1e1
#     assert rtf.size == 100
#     assert abs(rtf.r0 - 0.3) < 1e-10
#     assert rtf.alpha > 0
#     assert isinstance(rtf.get_default_int1d(), SimpsonIntegrator1D)


def test_power_properties():
    cases = get_power_cases()
    for rmin, rmax in cases:
        rtf = PowerRadialTransform(rmin, rmax, 100)
        assert rtf.rmin == rmin
        assert rtf.rmax == rmax
        assert rtf.size == 100
        assert rtf.power >= 2
        assert isinstance(rtf.get_default_int1d(), StubIntegrator1D)


def test_identiy_string():
    rtf1 = IndentityRadialTransform(45)
    s = rtf1.to_string()


# def test_linear_string():
#     rtf1 = LinearRTransform(np.random.uniform(
#         1e-5, 5e-5), np.random.uniform(1, 5), 88)
#     s = rtf1.to_string()
#     rtf2 = RTransform.from_string(s)
#     assert rtf1.rmin == rtf2.rmin
#     assert rtf1.rmax == rtf2.rmax
#     assert rtf1.size == rtf2.size
#     assert rtf1.alpha == rtf2.alpha

#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('LinearRTransform A 5')
#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('LinearRTransform A 5 .1')

#     rtf3 = RTransform.from_string('LinearRTransform -1.0 12.15643216847 77')
#     assert rtf3.rmin == -1.0
#     assert rtf3.rmax == 12.15643216847
#     assert rtf3.size == 77
#     assert rtf3.alpha > 0


# def test_exp_string():
#     rtf1 = ExpRTransform(np.random.uniform(1e-5, 5e-5),
#                          np.random.uniform(1, 5), 111)
#     s = rtf1.to_string()
#     rtf2 = RTransform.from_string(s)
#     assert rtf1.rmin == rtf2.rmin
#     assert rtf1.rmax == rtf2.rmax
#     assert rtf1.size == rtf2.size
#     assert rtf1.alpha == rtf2.alpha

#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('ExpRTransform A 5')
#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('ExpRTransform A 5 .1')

#     rtf3 = RTransform.from_string('ExpRTransform 1.0 12.15643216847 5')
#     assert rtf3.rmin == 1.0
#     assert rtf3.rmax == 12.15643216847
#     assert rtf3.size == 5
#     assert rtf3.alpha > 0


# def test_shifted_exp_string():
#     rtf1 = ShiftedExpRTransform(np.random.uniform(
#         1e-4, 5e-4), np.random.uniform(1e-5, 5e-5), np.random.uniform(1, 5), 781)
#     s = rtf1.to_string()
#     rtf2 = RTransform.from_string(s)
#     assert rtf1.rmin == rtf2.rmin
#     assert rtf1.rshift == rtf2.rshift
#     assert rtf1.rmax == rtf2.rmax
#     assert rtf1.size == rtf2.size
#     assert rtf1.r0 == rtf2.r0
#     assert rtf1.alpha == rtf2.alpha

#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('ShiftedExpRTransform A 5')
#     with assert_raises(ValueError):
#         rtf3 = RTransform.from_string('ShiftedExpRTransform A 5 .1 14')

#     rtf3 = RTransform.from_string(
#         'ShiftedExpRTransform 1.0 0.5 12.15643216847 5')
#     assert rtf3.rmin == 1.0
#     assert rtf3.rshift == 0.5
#     assert rtf3.rmax == 12.15643216847
#     assert rtf3.size == 5
#     assert rtf3.alpha > 0


def test_power_string():
    rtf1 = PowerRadialTransform(np.random.uniform(1e-3, 5e-3),
                                np.random.uniform(2.0, 5.0), 11)
    s = rtf1.to_string()
