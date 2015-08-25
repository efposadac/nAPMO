/*file: cuda_helper.cuh
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

/*Helper function to check errors in each CUDA function excecution.*/
#define CUERR                                                                              \
  {                                                                                        \
    cudaError_t err;                                                                       \
    if ((err = cudaGetLastError()) != cudaSuccess)                                         \
    {                                                                                      \
      printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
      printf("Aborting...\n");                                                      \
      exit(-1);                                                                            \
    }                                                                                      \
  }

#ifdef __CUDACC__
/*
Implementation of atomicAdd for double precision variables.
*/
__device__ __forceinline__ double atomicAdd(double *address, double val)
{
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do
  {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

/*
Implementation of double2 multiply operation.
*/
__device__ __forceinline__ double2 mult_double2(const double2 &a, const double2 &b)
{
  double2 r;
  r.x = a.x * b.x;
  r.y = a.y * b.y;
  return r;
}

#endif

#endif

