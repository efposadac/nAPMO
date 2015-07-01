/*file: gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "gauss_chebyshev.h"

void gaussChebyshev(int n, double* abscissas, double* weights) {
  // #ifdef _CUDA
  // 	gaussChebyshev_cuda(n, abscissas, weights);
  // #else
  int i;

  double aux = M_PI * 1.0 / (n + 1.0);

#ifdef _OMP
#pragma omp parallel for
#endif
  for (i = 0; i < n; ++i) {
    abscissas[i] = cos((i + 1) * aux);
    weights[i] = (aux * (1.0 - abscissas[i] * abscissas[i]));
  }
  // #endif
}