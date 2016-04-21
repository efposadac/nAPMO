/*file: gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "include/gauss_chebyshev.h"

void gaussChebyshev(int n, double rm, double *abscissas, double *weights) {
  int i;
  double aux_a, aux_b, aux_w;
  double aux = M_PI / (n + 1.0);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, aux_a, aux_w)
#endif
  for (i = 0; i < n; ++i) {
    aux_a = cos((i + 1) * aux);
    aux_b = sin((i + 1) * aux);

    aux_w = aux * aux_b * aux_b;

    /* Scale from interval (-1, 1) to (0, inf) */
    abscissas[i] = rm * (1.0 + aux_a) / (1.0 - aux_a); // (25)
    // divide sqrt(1-x^2) because Chebyshev second kind.
    weights[i] = aux_w / sqrt(1.0 - aux_a * aux_a);
    weights[i] *= 2.0 / ((aux_a - 1.0) * (aux_a - 1.0)); // (25)'
  }
}
