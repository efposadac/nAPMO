/*file: atomic_grid.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

// TODO: Update with the changes on the serial code....

extern "C" {
#include "include/atomic_grid.h"
}

#define THREADS_PER_BLOCK 64

double atomic_grid_integrate(AtomicGrid *grid, const int segments, double *f) {
  double *work, *weights_d, *f_d;
  double integral, *integral_d;

  int size = grid->size;

  cudaMalloc((void **)&work, size * sizeof(double));

  if (segments > 1) {
    cudaMalloc((void **)&f_d, size * segments * sizeof(double));
    cudaMemcpy(f_d, f, size * segments * sizeof(double),
               cudaMemcpyHostToDevice);
    dim3 dimGrid(((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1, 1);
    multiply_segmented_array<<<dimGrid, THREADS_PER_BLOCK>>>(size, segments,
                                                             f_d, work);
    CUERR
    cudaFree(f_d);
  } else {
    cudaMemcpy(work, f, size * sizeof(double), cudaMemcpyHostToDevice);
  }

  /*Initialize integral Value to 0.0*/
  cudaMalloc((void **)&integral_d, sizeof(double));
  cudaMemset(integral_d, 0.0, sizeof(double));

  cudaMalloc((void **)&weights_d, size * sizeof(double));
  cudaMemcpy(weights_d, grid->weights, size * sizeof(double),
             cudaMemcpyHostToDevice);

  dim3 dimGrid(((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1, 1);
  atomic_grid_integrate_kernel<<<dimGrid, THREADS_PER_BLOCK>>>(
      size, segments, work, weights_d, integral_d);

  integral = 0.0;
  cudaMemcpy(&integral, integral_d, sizeof(double), cudaMemcpyDeviceToHost);

  cudaFree(work);
  cudaFree(weights_d);
  cudaFree(integral_d);

  return integral * 4.0 * M_PI * grid->radii;
}

__global__ void atomic_grid_integrate_kernel(const int size, const int segments,
                                             double *work, double *weights,
                                             double *integral) {

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  __shared__ double temp[THREADS_PER_BLOCK];
  __shared__ double sum_block;

  temp[threadIdx.x] = 0.0;
  sum_block = 0.0;

  if (i < size) {
    temp[threadIdx.x] += work[i] * weights[i];
  }
  __syncthreads();
  atomicAdd(&sum_block, temp[threadIdx.x]);

  __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(integral, sum_block);
  }
}

__global__ void multiply_segmented_array(const int size, const int segments,
                                         double *f, double *output) {

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  __shared__ double temp[THREADS_PER_BLOCK];

  int j;
  if (i < size) {
    temp[threadIdx.x] = f[i];
    for (j = 1; j < segments; ++j) {
      temp[threadIdx.x] *= f[j * size + i];
    }
    output[i] = temp[threadIdx.x];
  }
}