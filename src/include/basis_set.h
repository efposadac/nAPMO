/*file: basis.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef BASIS_SET_H
#define BASIS_SET_H

struct _basis_set {
  int n_cont;            // Number of contractions.
  int *n_prim_cont;      // Number of primitives for each contraction.
  int *basis_l;          // Angular moment index sizeof (3*n_cont).
  double *exponent;      // exponent of each primitive.
  double *coefficient;   // contraction coefficient of each primitive.
  double *normalization; // Normalization constant for each primitive.
  double *origin;        // origin for each contraction.
};

typedef struct _basis_set BasisSet;

#ifdef _CUDA
#include "cuda_helper.cuh"
/*
Copy the host BasisSet structure into the device.
*/
void basis_set_init(BasisSet *basis, BasisSet *basis_d);

/*
Free the memory used on the CUDA device.
*/
void basis_set_free(BasisSet *basis_d);

#ifdef __CUDACC__
/*
Calculate the basis-set at point r
*/
/*
Calculate the basis-set at point r
*/
__device__ __forceinline__ void basis_set_compute_gto(BasisSet basis, double *r,
                                      double *output) {
  const int n_cont = basis.n_cont;
  double factor, RP2, temp, function_value;

  int counter = 0;
  for (int i = 0; i < n_cont; ++i) {

    int aux = i * 3;
    factor = 1.0, RP2 = 0.0;

    for (int j = 0; j < 3; ++j) {
      temp = r[j] - basis.origin[aux + j];
      RP2 += (temp * temp);
      factor *= pow(temp, basis.basis_l[aux + j]);
    }

    function_value = 0.0;
    for (int j = 0; j < basis.n_prim_cont[i]; ++j) {
      function_value +=
          basis.coefficient[counter] * exp(-basis.exponent[counter] * RP2);
      counter += 1;
    }

    function_value *= factor * basis.normalization[i];
    output[i] = function_value;
  }
}

#endif

#endif

#endif
