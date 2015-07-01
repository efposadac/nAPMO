/*file: gauss_chebyshev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef GAUSS_CHEBYSHEV_CUH
#define GAUSS_CHEBYSHEV_CUH

#include <math.h>
#include <stdio.h>

#define THREADS_PER_BLOCK 64

#define double_t float
/*
Computes abscissas and weights for the Gauss-Chebyshev quadrature of second
kind.
*/
void gaussChebyshev_cuda(int n, double_t *abscissas, double_t *weights);

#endif