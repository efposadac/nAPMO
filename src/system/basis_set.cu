/*file: basis_set.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

#define THREADS_PER_BLOCK 64

#include "basis_set.h"

#include "cuda_helper.cuh"

void basis_set_init(BasisSet *basis, BasisSet *basis_d) {
  int i, basisSize;

  /*Calculating total size of the basis set*/
  basisSize = 0;
  for (i = 0; i < basis->n_cont; ++i) {
    basisSize += basis->n_prim_cont[i];
  }

  int bytes_int = basis->n_cont * sizeof(int);
  int bytes_double = basis->n_cont * sizeof(double);
  int bytes_basis = basisSize * sizeof(double);

  /*Allocating space for the device structure*/
  cudaMalloc((void **)&basis_d->n_prim_cont, bytes_int);
  cudaMalloc((void **)&basis_d->prim_index, bytes_int);
  cudaMalloc((void **)&basis_d->basis_l, 3 * bytes_int);
  cudaMalloc((void **)&basis_d->origin, 3 * bytes_double);
  cudaMalloc((void **)&basis_d->normalization, bytes_double);
  cudaMalloc((void **)&basis_d->exponent, bytes_basis);
  cudaMalloc((void **)&basis_d->coefficient, bytes_basis);

  /*Copying data to device*/
  basis_d->n_cont = basis->n_cont;
  cudaMemcpy(basis_d->n_prim_cont, basis->n_prim_cont, bytes_int,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->prim_index, basis->prim_index, bytes_int,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->basis_l, basis->basis_l, 3 * bytes_int,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->normalization, basis->normalization, bytes_double,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->origin, basis->origin, 3 * bytes_double,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->exponent, basis->exponent, bytes_basis,
             cudaMemcpyHostToDevice);
  cudaMemcpy(basis_d->coefficient, basis->coefficient, bytes_basis,
             cudaMemcpyHostToDevice);
}

void basis_set_free(BasisSet *basis_d) {
  cudaFree(basis_d->n_prim_cont);
  cudaFree(basis_d->prim_index);
  cudaFree(basis_d->basis_l);
  cudaFree(basis_d->normalization);
  cudaFree(basis_d->origin);
  cudaFree(basis_d->exponent);
  cudaFree(basis_d->coefficient);
}
