/*file: density.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/density.h"

#ifndef _CUDA

void density_gto(BasisSet *basis, double *r, double *dens, double *output,
                 int size) {
  int point, idx;

#ifdef _OMP
#pragma omp parallel for default(shared) private(point, idx)
#endif
  for (point = 0; point < size; ++point) {
    idx = point * 3;
    output[point] = density_gto_r(basis, &r[idx], dens);
  }
}

double density_gto_r(BasisSet *basis, double r[3], double *dens) {
  int i, j, idx, counter = 0;
  int n_cont = basis->n_cont;
  double temp;
  double function_value, output = 0.0;
  double RP2, factor;

  double *basis_val = (double *)calloc(n_cont * 2, sizeof(double));

  for (i = 0; i < n_cont; ++i) {
    factor = 1.0, RP2 = 0.0;
    for (j = 0; j < 3; ++j) {
      idx = i * 3;
      temp = r[j] - basis->origin[idx + j];
      RP2 += (temp * temp);
      factor *= pow(temp, basis->basis_l[idx + j]);
    }

    function_value = 0.0;
    for (j = 0; j < basis->n_prim_cont[i]; ++j) {
      function_value +=
          basis->coefficient[counter] * exp(-basis->exponent[counter] * RP2);
      counter += 1;
    }

    function_value *= factor * basis->normalization[i];

    for (j = 0; j < n_cont; ++j) {
      basis_val[n_cont + j] += function_value * dens[i * n_cont + j];
    }

    basis_val[i] = function_value;
  }

  for (i = 0; i < n_cont; ++i) {
    output += basis_val[i] * basis_val[n_cont + i];
  }

  free(basis_val);
  return output;
}
#endif