# test.py

import numpy as np
from napmo.utilities.radial_quadratures import *

# 3x^2+2


def print_matrix(A, n):
    for i in range(n):
        for j in range(n):
            print("%12.5f" % (A[i, j]), end="")
        print("")


def f(x):
    return 3*x**2+2


def Df(x):
    return 6*x


# Ax=b
def first_2p(n, b):
    np2 = n + 2

    x_i = chebgauss_z(np2)
    h = x_i[0]

    A = np.zeros([np2, np2])
    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    for i in range(1, n+1):
        A[i, i-1] = -0.5/h
        A[i, i+1] = 0.5/h

    x = np.linalg.solve(A, Df_a)

    return x


def first_4p_f(f, i, x, h):
    aux = (
        - (25.0 / 12.0) * f(x[i]) +
        4.0 * f(x[i+1]) -
        3.0 * f(x[i+2]) +
        (4.0 / 3.0) * f(x[i+3]) -
        (1.0 / 4.0) * f(x[i+4])
    )

    return aux / h


def first_4p_b(f, i, x, h):
    aux = (
        (25.0 / 12.0) * f(x[i]) -
        4.0 * f(x[i-1]) +
        3.0 * f(x[i-2]) -
        (4.0 / 3.0) * f(x[i-3]) +
        (1.0 / 4.0) * f(x[i-4])
    )

    return aux / h


def first_7p_Bickley(f, i, x, h):
    aux = 1.0/(60.0*h)
    d = (
        f(x[i+3]) -
        9.0 * f(x[i+2]) +
        45.0 * f(x[i+1]) -
        45.0 * f(x[i-1]) +
        9.0 * f(x[i-2]) -
        f(x[i-3])
    )

    return d * aux


def first_7p(n):
    np6 = n + 6
    x_i = chebgauss_z(n+6)
    b = Df(x_i)

    # Boundary conditions
    b[0] = f(x_i[0])
    # b[1] = f(x_i[1])
    # b[2] = f(x_i[2])

    b[-1] = f(x_i[-1])
    # b[-2] = f(x_i[-2])
    # b[-3] = f(x_i[-3])

    h = x_i[0]

    A = np.zeros([np6, np6])

    A[0, 0] = 1.0
    A[-1, -1] = 1.0

    aux = 1.0 / (12.0*h)

    A[1, 0] = -25 * aux
    A[1, 1] = 48 * aux
    A[1, 2] = -36 * aux
    A[1, 3] = 16 * aux
    A[1, 4] = -3 * aux

    A[-2, -1] = 25 * aux
    A[-2, -2] = -48 * aux
    A[-2, -3] = 36 * aux
    A[-2, -4] = -16 * aux
    A[-2, -5] = 3 * aux

    aux = 1.0 / (60.0*h)

    A[2, 0] = -137 * aux
    A[2, 1] = 300 * aux
    A[2, 2] = -300 * aux
    A[2, 3] = 200 * aux
    A[2, 4] = -75 * aux
    A[2, 5] = 12 * aux

    A[-3, -1] = 137 * aux
    A[-3, -2] = -300 * aux
    A[-3, -3] = 300 * aux
    A[-3, -4] = -200 * aux
    A[-3, -5] = 75 * aux
    A[-3, -6] = -12 * aux

    for i in range(3, n+3):
        A[i, i-3] = -aux
        A[i, i-2] = 9.0 * aux
        A[i, i-1] = -45.0 * aux
        A[i, i+1] = 45.0 * aux
        A[i, i+2] = -9.0 * aux
        A[i, i+3] = aux

    print_matrix(A, np6)

    x = np.linalg.solve(A, b)

    return x, A, b


if __name__ == '__main__':
    n = 3

    f_n, A, b = first_7p(n)

    # Check output
    x_i = chebgauss_z(n+6)
    Df_a = Df(x_i)

    # print(x_i)

    f_a = f(x_i)

    print(np.abs(f_n - f_a))
    print(np.allclose(f_n, f_a))
    print(np.allclose(np.dot(A, f_n)[3:-3], Df_a[3:-3]))
    print("f_a", f_a)
    print("f_n", f_n)
    print("b", b)
    print("Dot", np.dot(A, f_n), Df_a)
