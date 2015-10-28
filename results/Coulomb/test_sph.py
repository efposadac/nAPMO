import numpy as np
import scipy as sci
import scipy.special as sp
from ctypes import *
import os

lp = os.path.join(os.path.dirname(__file__), 'test_gsl.so')
test_gsl = CDLL(lp)

test_gsl.real_sph.restype = c_double
test_gsl.real_sph.argtypes = [
    c_int,
    c_int,
    c_double,
    c_double
]


def real_spherical_harmonics(m, l, theta, phi):
    if m == 0:
        aux = sp.sph_harm(m, l, theta, phi)
    elif m < 0:
        aux_a = sp.sph_harm(-m, l, theta, phi)
        aux_b = sp.sph_harm(m, l, theta, phi)
        aux = 1.0j * 0.70710678118654757 * (aux_b + aux_a)
    else:
        aux_a = sp.sph_harm(m, l, theta, phi)
        aux_b = sp.sph_harm(-m, l, theta, phi)
        aux = 0.70710678118654757 * (aux_b - aux_a)
    return aux.real

lmax = 4
for l in range(lmax+1):
    for m in range(-l, l+1):
        print(real_spherical_harmonics(m, l, 0.0, 0.0), test_gsl.real_sph(l, m, 0.0, 0.0))
