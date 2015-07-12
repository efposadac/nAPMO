/*file: cuda_helper.cuh
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef CUDA_HELPER_H
#define CUDA_HELPER_H

#include <stdio.h>
#include <cuda_runtime.h>

/*Helper function to check errors in each CUDA function excecution.*/
#define CUERR                                                                              \
  {                                                                                        \
    cudaError_t err;                                                                       \
    if ((err = cudaGetLastError()) != cudaSuccess)                                         \
    {                                                                                      \
      printf("CUDA error: %s, %s line %d\n", cudaGetErrorString(err), __FILE__, __LINE__); \
      printf("Thread aborting...\n");                                                      \
    }                                                                                      \
  }

#endif