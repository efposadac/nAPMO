/*file: becke_grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef BECKE_GRID_H
#define BECKE_GRID_H

#include <stdlib.h>
#include <stdio.h>

#ifdef _OMP
#include <omp.h>
#endif

#include "lebedev.h"
#include "gauss_chebyshev.h"
#include "system.h"

struct _grid
{
    int n_radial;              // Number of radial points.
    int n_angular;             // Number of angular points.
    double* radial_abscissas;  // Radial abscissas (len = n_radial)
    double* radial_weights;    // Radial weights (len = n_radial)
    double* angular_theta;     // theta coordinate of angular quadrature (len =
                               // n_angular)
    double* angular_phi;       // phi coordinate of angular quadrature (len = n_angular)
    double* angular_weights;   // Weights of angular quadrature (len = n_angular)
};
typedef struct _grid Grid;

/*
Initialization of grid. Memory allocation in C.
Calculates the n_radial Gauss-Chebyshev quadrature and n_angular Lebedev
quadrature points.
References:
    Becke, A. D. A multicenter numerical integration scheme for polyatomic
molecules. J. Chem. Phys. 88, 2547 (1988).
*/
void grid_init(Grid* grid);

/*
Free memory used in grid construction.
*/
void grid_free(Grid* grid);

/*
Computes the Becke weights :math:`w(r)` at point ``r`` for particle
``particleID`` as described in eq. 22 Becke, 1988.

References:
    Becke, A. D. A multicenter numerical integration scheme for polyatomic
molecules. J. Chem. Phys. 88, 2547 (1988).

Args:
    (double[3]): Point of the grid in which the weight will be calculated.
    particleID (int): The particle index who owns the ``r`` point.
    sys (System): system structure.

Returns:
    output (double): The value of cell_function (eq. 13, Becke, 1988) at point
``r``
*/
double grid_weights(System* sys, double r[3], int particleID);

/*
Functional :math:`\rho({\\bf r})`
*/
double grid_integrate(System* sys, Grid* grid);

/*
Functional :math:`\rho({\\bf r})`
*/
double grid_density(System* sys, double* r, double* dens);

/*
Iterated cutoff profile. eq. 21, Becke 1988.
*/
static inline double grid_soft_mu(const double mu) { return 0.5 * mu * (3.0 - mu * mu); }

#endif