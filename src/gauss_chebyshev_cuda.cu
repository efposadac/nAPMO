/*file: gauss_chebyshev_cuda.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include <cuda_runtime.h>

extern "C"{
#include "gauss_chebyshev_cuda.cuh"
}

__global__ void kernel_gaussChebyshev_cuda(const int n, double_t *abscissas,
                                           double_t *weights) {
  const unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
  double_t aux_1 =  M_PI / (n + 1.0);
  double_t aux_2;

  if (i < n) {
    aux_2 = cos((i + 1) * aux_1);
    weights[i] = (aux_1 * (1.0 - aux_2 * aux_2));
    abscissas[i] = aux_2;
  }
}

void gaussChebyshev_cuda(int n, double_t *abscissas, double_t *weights) {
  double_t *g_a;
  double_t *g_w;

  g_a = g_w = NULL;

  cudaMalloc((void **)&g_a, n * sizeof(double_t));
  cudaMalloc((void **)&g_w, n * sizeof(double_t));

  int gridDim = ((n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);

  kernel_gaussChebyshev_cuda <<<gridDim, THREADS_PER_BLOCK>>> (n, g_a, g_w);

  /* check if the kernel launch was successful */
  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    printf("failed to lauch GPU kernel:\n%s\n", cudaGetErrorString(err));
    return;
  }

  cudaMemcpy(abscissas, g_a, n * sizeof(double_t), cudaMemcpyDeviceToHost);
  cudaMemcpy(weights, g_w, n * sizeof(double_t), cudaMemcpyDeviceToHost);

  cudaFree(g_a);
  cudaFree(g_w);
}
