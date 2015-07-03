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

#ifdef _CUDA
#include "gauss_chebyshev_cuda.cuh"
#endif

/*
Computes abscissas and weights for the Gauss-Chebyshev quadrature of second
kind.
*/
void gaussChebyshev(int n, double* abscissas, double* weights);


#endif