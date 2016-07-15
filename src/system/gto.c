/*file: gto.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "gto.h"

double gto_normalize_primitive(PrimitiveGaussian *f) {
  int i;
  double output;
  int aux = 0;

  for (i = 0; i < 3; ++i) {
    aux += f->l[i];
  }
  output = (pow((2.0 * f->exponent / M_PI), 0.75)) /
           sqrt(utils_factorial2(abs(2 * f->l[0] - 1)) *
                utils_factorial2(abs(2 * f->l[1] - 1)) *
                utils_factorial2(abs(2 * f->l[2] - 1)) /
                (pow((4.0 * f->exponent), aux)));

  return output;
}

void gto_compute_primitive(PrimitiveGaussian *f, double *coord, double *output,
                           const int n_coord) {
  int i, j, idx;
  double RP2 = 0.0, factor = 1.0, aux;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, j, idx, aux, factor, RP2)
#endif
  for (i = 0; i < n_coord; ++i) {
    RP2 = 0.0;
    factor = 1.0;
    idx = i * 3;
    for (j = 0; j < 3; ++j) {
      aux = coord[idx + j] - f->origin[j];
      factor *= pow(aux, f->l[j]);
      RP2 += aux * aux;
    }
    output[i] =
        f->coefficient * f->normalization * factor * exp(-f->exponent * RP2);
  }
}
