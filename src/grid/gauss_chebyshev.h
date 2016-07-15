/*file: gauss_chebyshev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/
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

#endif