/*file: gauss_chebyshev_cuda.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "gauss_chebyshev.h"
}

#include "cuda_helper.cuh"

#define THREADS_PER_BLOCK 64

__global__ void gaussChebyshev_cuda_kernel(const int n, double *abscissas, double *weights)
{
  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  
  const double aux_1 = M_PI / (n + 1.0);
  double aux_2;

  if (i < n)
  {
    aux_2 = cos((i + 1) * aux_1);
    weights[i] = (aux_1 * (1.0 - aux_2 * aux_2));
    abscissas[i] = aux_2;
  }
}

void gaussChebyshev_cuda(const int n, double *abscissas, double *weights)
{
  int gridDim = ((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  gaussChebyshev_cuda_kernel <<<gridDim, THREADS_PER_BLOCK>>> (n, abscissas, weights);

  /* check if the kernel launch was successful */
  CUERR
}
