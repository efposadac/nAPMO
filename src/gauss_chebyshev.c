/*file: gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "gauss_chebyshev.h"

void gaussChebyshev(int n, double rm, double *abscissas, double *weights) {
  int i;
  double aux_a, aux_w;
  double aux = M_PI / (n + 1.0);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, aux_a, aux_w)
#endif
  for (i = 0; i < n; ++i) {
    aux_a = cos((i + 1) * aux);
    aux_w = aux * (1.0 - aux_a * aux_a);

    /* Scale from interval (-1, 1) to (0, inf) */
    abscissas[i] = rm * (1.0 + aux_a) / (1.0 - aux_a); // (25)
    weights[i] =
        aux_w /
        sqrt(1.0 -
             aux_a *
                 aux_a); // divide sqrt(1-x^2) because Chebyshev second kind.
    weights[i] *= 2.0 / ((aux_a - 1.0) * (aux_a - 1.0)); // (25)'
  }
}

void gaussChebyshev_get_z(int n, double rm, double *abscissas, double *z) {

  int i;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < n; ++i) {
    z[i] = acos((abscissas[i] - rm) / (abscissas[i] + rm)) * M_1_PI;
  }
}

void gaussChebyshev_deriv_z(int n, double rm, double *abscissas,
                            double *deriv_z) {
  int i;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < n; ++i) {
    deriv_z[i] = -sqrt(rm * abscissas[i] / pow(rm + abscissas[i], 2)) /
                 (M_PI * abscissas[i]);
  }
}

void gaussChebyshev_deriv2_z(int n, double rm, double *abscissas,
                             double *deriv2_z) {
  int i;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < n; ++i) {
    deriv2_z[i] =
        rm * rm * (rm + (3.0 * abscissas[i])) /
        (2.0 * M_PI * pow(rm * abscissas[i] / pow(rm + abscissas[i], 2), 1.5) *
         pow(rm + abscissas[i], 5));
  }
}
