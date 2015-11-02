/*file: gauss_chebyshev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef GAUSS_CHEBYSHEV_H
#define GAUSS_CHEBYSHEV_H

#include <math.h>

#ifdef _OMP
#include <omp.h>
#endif

/*
Computes abscissas and weights for the Gauss-Chebyshev quadrature of second
kind. Includes the interval transformation from [-1,1] to [0,inf] according to
Becke's paper.
*/
void gaussChebyshev(int n, double rm, double *abscissas, double *weights);

/*
Returns the radial points mapped uniformly in the interval [0,1], see
Becke's paper.
*/
void gaussChebyshev_get_z(int n, double rm, double *abscissas, double *z);

/*
Returns the first derivative of the uniform z grid.
*/
void gaussChebyshev_deriv_z(int n, double rm, double *abscissas,
                            double *deriv_z);

/*
Returns the second derivative of the uniform z grid.
*/
void gaussChebyshev_deriv2_z(int n, double rm, double *abscissas,
                             double *deriv2_z);

#endif