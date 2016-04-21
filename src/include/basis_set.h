/*file: basis.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/
#ifndef BASIS_SET_H
#define BASIS_SET_H

struct _basis_set {
  int n_cont;            // Number of contractions.
  int *n_prim_cont;      // Number of primitives for each contraction.
  int *prim_index;       // Indexes to access to primitives arrays.
  int *basis_l;          // Angular moment index, size (3*n_cont).
  double *origin;        // origin for each contraction.
  double *normalization; // Normalization constant for each primitive.
  double *exponent;      // exponent of each primitive.
  double *coefficient;   // contraction coefficient of each primitive.
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
Calculate the contraction ``cont`` at point ``r``
*/
__device__ __forceinline__ double basis_set_compute_gto(BasisSet basis,
                                                        int cont, double *r) {

  const int n_prim = basis.n_prim_cont[cont];
  const int prim_index = basis.prim_index[cont];
  const int aux = cont * 3;

  double factor, RP2, temp, function_value;

  factor = 1.0, RP2 = 0.0;
  for (int j = 0; j < 3; ++j) {
    temp = r[j] - basis.origin[aux + j];
    RP2 += (temp * temp);
    factor *= pow(temp, basis.basis_l[aux + j]);
  }

  function_value = 0.0;
  for (int j = 0; j < n_prim; ++j) {
    function_value += basis.coefficient[prim_index + j] *
                      exp(-basis.exponent[prim_index + j] * RP2);
  }

  function_value *= factor * basis.normalization[cont];

  return function_value;
}

#endif

#endif

#endif
