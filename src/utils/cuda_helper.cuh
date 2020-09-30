/*file: cuda_helper.cuh
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

/*Helper function to check errors in each CUDA function excecution.*/
#define CUERR                                                                  \
  {                                                                            \
    cudaError_t err;                                                           \
    if ((err = cudaGetLastError()) != cudaSuccess) {                           \
      printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err),          \
             __FILE__, __LINE__);                                              \
      printf("Aborting...\n");                                                 \
      exit(-1);                                                                \
    }                                                                          \
  }

#ifdef __CUDACC__
/*
Implementation of atomicAdd for double precision variables.
*/
__device__ __forceinline__ double atomicAdd(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

/*
Implementation of double2 multiply operation.
*/
__device__ __forceinline__ double2 mult_double2(const double2 &a,
                                                const double2 &b) {
  double2 r;
  r.x = a.x * b.x;
  r.y = a.y * b.y;
  return r;
}

__device__ __forceinline__ double distance(double a[3], double b[3]) {
  int i;
  double output, aux;

  output = 0.0;
  for (i = 0; i < 3; ++i) {
    aux = a[i] - b[i];
    output += aux * aux;
  }

  return sqrt(output);
}

#endif

#endif
