# file: test_extrapolation.py
# nAPMO package
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

from napmo.utilities.cubic_spline import CubicSpline
from napmo.utilities.extrapolation import *

import numpy as np


def test_extrapolation1_identity():
    x = np.arange(10, dtype=float)
    y = np.exp(-0.3 * x)
    d = -0.3 * y
    cs = CubicSpline(y, d)
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx) - np.exp(-0.3 * newx)).max() < 1e-10
    assert abs(cs.deriv(newx) - -0.3 * np.exp(-0.3 * newx)).max() < 1e-10
    newx = np.array([10.5, 11.5])
    assert abs(cs(newx)).max() < 1e-10
    assert abs(cs.deriv(newx)).max() < 1e-10


def test_extrapolation2_identity():
    x = np.arange(10, dtype=float)
    y = x**2 + 1
    d = 2 * x
    cs = CubicSpline(y, d)
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx) - 1.0).max() < 1e-10
    assert abs(cs.deriv(newx)).max() < 1e-10


# def test_extrapolation1_exp():
#     rtf = ExpRTransform(0.1, 1.0, 10)
#     x = rtf.get_radii()
#     y = np.exp(-0.3 * x)
#     d = -0.3 * y
#     cs = CubicSpline(y, d, rtf)
#     newx = np.array([0.001, 0.01])
#     assert abs(cs(newx) - np.exp(-0.3 * newx)).max() < 1e-10
#     assert abs(cs.extrapolation.eval_left(
#         newx[0]) - np.exp(-0.3 * newx[0])).max() < 1e-10
#     assert abs(cs.deriv(newx) - -0.3 * np.exp(-0.3 * newx)).max() < 1e-10
#     assert abs(cs.extrapolation.deriv_left(
#         newx[0]) - -0.3 * np.exp(-0.3 * newx[0])).max() < 1e-10
#     newx = np.array([1.1, 10.1])
#     assert abs(cs(newx)).max() < 1e-10
#     assert abs(cs.deriv(newx)).max() < 1e-10


# def test_bug_steven():
#     y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
#                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.95307674e-005,
#                   2.64520014e-004, 4.86536423e-004, 7.72847240e-004, 1.15087885e-003, 1.66508373e-003, 2.38709713e-003,
#                   3.42860950e-003, 4.95314131e-003, 7.17684800e-003, 1.03407874e-002, 1.46341247e-002, 2.00627694e-002,
#                   2.63026157e-002, 3.26349709e-002, 3.80678816e-002, 4.16324963e-002, 4.26735205e-002, 4.09493929e-002,
#                   3.65937608e-002, 3.01384965e-002, 2.25802551e-002, 1.52220562e-002, 9.20415761e-003, 5.04355509e-003,
#                   2.56723646e-003, 1.24504774e-003, 5.72945090e-004, 2.36808736e-004, 8.00848384e-005, 1.99672917e-005,
#                   3.30308018e-006, 3.22840493e-007, 1.62065171e-008, 3.50644056e-010, 2.62139812e-012, 5.12409163e-015,
#                   1.84238101e-018, 7.81842741e-023, 2.23721221e-028, 2.12999970e-035, 2.76742060e-044, 1.59368870e-055,
#                   9.84245631e-070, 1.08780738e-087, 2.24648735e-110])

#     rtf = ExpRTransform(0.0003779452267842504, 37.79452267842504, 100)
#     cs = CubicSpline(y, rtransform=rtf)
#     d = np.array([0.000141, 0.00141, 0.0141, 0.141, 1.41])
#     s = cs(d)
#     assert not np.isnan(s).any()


def test_extrapolation_identity_power():
    x = np.arange(10, dtype=float)
    y = 1 / (np.exp(-2 * x) + x**2)
    d = -y**2 * (np.exp(-2 * x) * (-2) + 2 * x)
    cs = CubicSpline(y, d, extrapolation=PowerExtrapolation(-2.0))
    assert cs.extrapolation.power == -2.0
    newx = np.array([-2.5, -1.1])
    assert abs(cs(newx)).max() == 0.0
    assert abs(cs.deriv(newx)).max() == 0.0
    newx = np.array([10.5, 11.5])
    assert abs(cs(newx) - 1 / newx**2).max() < 1e-10
    assert abs(cs.extrapolation.eval_right(
        newx[0]) - 1 / newx[0]**2).max() < 1e-10
    assert abs(cs.deriv(newx) - -2 / newx**3).max() < 1e-10
    assert abs(cs.extrapolation.deriv_right(newx[0]) - -2 / newx[0]**3) < 1e-10


# def test_extrapolation_exp_power():
#     rtf = ExpRTransform(0.1, 1.0, 10)
#     x = rtf.get_radii()
#     y = x**2
#     d = 2 * x
#     cs = CubicSpline(y, d, rtf, PowerExtrapolation(2))
#     assert cs.extrapolation.power == 2.0
#     newx = np.array([0.001, 0.01])
#     assert abs(cs(newx)).max() == 0.0
#     assert abs(cs.deriv(newx)).max() == 0.0
#     newx = np.array([1.1, 10.1])
#     assert abs(cs(newx) - newx**2).max() < 1e-10
#     assert abs(cs.extrapolation.eval_right(newx[0]) - newx[0]**2).max() < 1e-10
#     assert abs(cs.deriv(newx) - 2 * newx).max() < 1e-10
#     assert abs(cs.extrapolation.deriv_right(newx[0]) - 2 * newx[0]) < 1e-10


# def test_consistency_h5():
#     with h5.File('horton.grid.test.test_cubic_spline.test_consistency_h5', driver='core', backing_store=False) as chk:
#         rtf = ExpRTransform(0.1, 1.0, 10)
#         y = np.random.normal(0, 1, 10)
#         d = np.random.normal(0, 1, 10)
#         cs1 = CubicSpline(y, d, rtf, PowerExtrapolation(2))
#         cs1.to_hdf5(chk)
#         cs2 = CubicSpline.from_hdf5(chk)
#         assert (cs1.y == cs2.y).all()
#         assert (cs1.dx == cs2.dx).all()
#         assert (cs1.dt == cs2.dt).all()
#         assert cs1.rtransform.to_string() == cs2.rtransform.to_string()
#         assert cs1.extrapolation.to_string() == cs2.extrapolation.to_string()


# def test_zero_extrapolation_string():
#     assert ZeroExtrapolation().to_string() == 'ZeroExtrapolation'
#     assert isinstance(Extrapolation.from_string(
#         'ZeroExtrapolation'), ZeroExtrapolation)


# def test_cusp_extrapolation_string():
#     assert CuspExtrapolation().to_string() == 'CuspExtrapolation'
#     assert isinstance(Extrapolation.from_string(
#         'CuspExtrapolation'), CuspExtrapolation)


# def test_power_extrapolation_string():
#     assert PowerExtrapolation(2.0).to_string() == 'PowerExtrapolation 2.0'
#     ep = Extrapolation.from_string(
#         PowerExtrapolation(5.1247953315476).to_string())
#     assert isinstance(ep, PowerExtrapolation)
#     assert ep.power == 5.1247953315476
