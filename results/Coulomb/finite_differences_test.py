# test.py

import numpy as np
from napmo.utilities.radial_quadratures import *
import matplotlib.pyplot as plt

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
    b[-1] = f(x_i[-1])

    h = x_i[0]

    A = np.zeros([np6, np6])

    for i in range(1):
        A[i, i] = 1.0
        A[-i-1, -i-1] = 1.0

    aux = 1.0 / (2.0 * h)

    A[1, 0] = -aux
    A[1, 2] = aux

    A[-2, -3] = -aux
    A[-2, -1] = aux

    aux = 1.0 / (12.0*h)
    p5 = np.array([1.0, -8.0])
    for j in range(2):
        A[2, 2-2+j] = p5[j] * aux
        A[2, 2+2-j] = -p5[j] * aux
        A[-3, -3-2+j] = p5[j] * aux
        A[-3, -3+2-j] = -p5[j] * aux

    aux = 1.0 / (60.0*h)
    p7 = np.array([-1.0, 9.0, -45.0])
    for i in range(3, n+3):
        for j in range(3):
            A[i, i-3+j] = p7[j] * aux
            A[i, i+3-j] = -p7[j] * aux

    # print_matrix(A, np6)

    x = np.linalg.solve(A, b)

    return x, A, b


if __name__ == '__main__':
    n = 8

    f_n, A, b = first_7p(n)

    # Check output
    x_i = chebgauss_z(n+6)
    Df_a = Df(x_i)

    f_a = f(x_i)

    print(np.allclose(f_n, f_a))
    print(np.allclose(np.dot(A, f_n)[-1:1], Df_a[-1:1]))

    plt.plot(x_i, f_a, label="f Analytic")
    plt.plot(x_i, f_n, 'o', label="f Numeric")
    plt.plot(x_i, Df_a, label="D Analytic")
    plt.plot(x_i, A.dot(f_n), 'o', label="D Numeric")
    plt.legend(loc=2)
    plt.show()
