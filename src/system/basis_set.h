/*file: basis.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef BASIS_SET_H
#define BASIS_SET_H

#include <algorithm>
#include <math.h>

#include "gto.h"

struct BasisSet {
private:
  int nbasis;
  int max_l;
  int max_nprim;
  std::vector<ContractedGaussian> cont;

public:
  BasisSet() : nbasis(0), max_l(0), max_nprim(0){};

  BasisSet(const BasisSet &) = default;

  BasisSet(ContractedGaussian **contractions, int n);

  ~BasisSet(){};

  std::vector<double> compute(double *r);
  
  std::vector<double> deriv(double *r);

  void update(BasisSet *other);

  int get_nbasis() { return nbasis; };
  int get_max_l() { return max_l; };
  int get_max_nprim() { return max_nprim; };

  std::vector<ContractedGaussian> &get_cont() { return cont; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
BasisSet wrapper to python
*/

BasisSet *BasisSet_new_empty();

BasisSet *BasisSet_new(ContractedGaussian **primitives, int n);

void BasisSet_compute(BasisSet *basis, double *r, double *output, int size);

void BasisSet_deriv(BasisSet *basis, double *r, double *output, int size);

void BasisSet_update(BasisSet *basis, BasisSet *other);

int BasisSet_get_nbasis(BasisSet *basis);

int BasisSet_get_max_l(BasisSet *basis);

int BasisSet_get_max_nprim(BasisSet *basis);

#ifdef __cplusplus
}
#endif

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
                      basis.p_normalization[prim_index + j] *
                      exp(-basis.exponent[prim_index + j] * RP2);
  }

  function_value *= factor * basis.normalization[cont];

  return function_value;
}

#endif

#endif

#endif
