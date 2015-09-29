import numpy as np
import scipy as sci
import scipy.special as sp
import matplotlib.pyplot as plt

from napmo.utilities.angular_quadratures import *
from napmo.utilities.radial_quadratures import *
from napmo.interfaces.molecular_system import *
from napmo.interfaces.becke_grid import *


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


def Y(l, m, theta, phi):
    return real_spherical_harmonics(m, l, phi, theta)


# Define integrand P(r, theta, phi)
def int_l(theta, phi, func, r, l, m):

    # Convert to cartesian
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    coord = np.dstack((x, y, z)).reshape(len(x), 3)

    return func(coord) * Y(l, m, theta, phi)


def rho_lm(func, rad, n, lmax):
    # Size of expansion
    lsize = (lmax + 1) * (lmax + 1)

    # Start calculation
    p_lm = np.zeros((len(rad), lsize))

    for r in range(len(rad)):
        lindex = 0
        for l in range(lmax+1):
            for m in range(-l, l+1):
                # integrate
                p_lm[r, lindex] = lebedev_integrate(int_l, n, args=(func, rad[r], l, m))
                lindex += 1

    return p_lm


def recover_rho(grid, p_lm, lmax):
    idx = 0
    rho = np.zeros(grid.size)
    for r in range(grid.n_radial):
        for a in range(grid.n_angular):
            lindex = 0
            for l in range(lmax+1):
                for m in range(-l, l+1):
                    rho[idx] += p_lm[r, lindex] * Y(l, m, grid.angular_theta[a], grid.angular_phi[a])
                    lindex += 1
            idx += 1

    return rho

if __name__ == '__main__':

    # Molecule definition
    basis_file = "TEST.json"

    molecule = MolecularSystem()
    molecule.add_atom("He", [0.0, 0.0, 0.3704240745], basis_kind="GTO", basis_file=basis_file)
    molecule.show()

    particles = molecule.get('atoms')
    rm = particles[-1].get('atomic_radii_2')

    # Grid definition
    angularPoints = 110
    radialPoints = 50

    grid = BeckeGrid(radialPoints, angularPoints)
    grid.move(scaling_factor=rm)
    grid.show()

    # Functions
    basis = molecule.get_basis_set('e-')
    g_a = basis.get('function')[-1]

    # Problem
    lmax = int(lebedev_get_order(angularPoints)/2)
    rad_quad = np.array([grid.radial_abscissas[i] for i in range(grid.n_radial)]) * rm

    #########################################################
    def p_ab(r):
        return g_a.compute(r) * g_a.compute(r)

    p_rtp = p_ab(grid.xyz)

    #########################################################

    expansion = rho_lm(p_ab, rad_quad, angularPoints, lmax)

    #########################################################
    recovered = recover_rho(grid, expansion, lmax)

    plt.plot(rad_quad, p_rtp[::angularPoints], 'r-', label='Patron')
    plt.plot(rad_quad, recovered[::angularPoints], 'go', label='Recovered')
    plt.legend()
    plt.xlim([0, 3])
    plt.show()
    grid.free()
