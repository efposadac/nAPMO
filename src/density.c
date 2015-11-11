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
  int i, j, idx, aux, counter, point;
  int n_cont = basis->n_cont;
  double temp, temp_val = 0.0;
  double function_value;
  double RP2, factor;

#ifdef _OMP
#pragma omp parallel default(shared) private(                                  \
    point, idx, counter, j, aux, factor, RP2, temp, function_value)            \
        reduction(+ : temp_val)
#endif
  {
    double *basis_val = (double *)calloc(n_cont * 2, sizeof(double));
#ifdef _OMP
#pragma omp for
#endif
    for (point = 0; point < size; ++point) {
      idx = point * 3;
      counter = 0;
      temp_val = 0.0;

      for (j = 0; j < n_cont * 2; ++j) {
        basis_val[j] = 0.0;
      }

      for (i = 0; i < n_cont; ++i) {
        aux = i * 3;
        factor = 1.0; 
        RP2 = 0.0;
        for (j = 0; j < 3; ++j) {
          temp = r[idx + j] - basis->origin[aux + j];
          RP2 += (temp * temp);
          factor *= pow(temp, basis->basis_l[aux + j]);
        }

        function_value = 0.0;
        for (j = 0; j < basis->n_prim_cont[i]; ++j) {
          function_value += basis->coefficient[counter] *
                            exp(-basis->exponent[counter] * RP2);
          counter += 1;
        }

        function_value *= factor * basis->normalization[i];

        for (j = 0; j < n_cont; ++j) {
          basis_val[n_cont + j] += function_value * dens[i * n_cont + j];
        }

        basis_val[i] = function_value;
      }

      for (i = 0; i < n_cont; ++i) {
        temp_val += basis_val[i] * basis_val[n_cont + i];
      }
      output[point] = temp_val;
    }

    free(basis_val);
  }
}

#endif