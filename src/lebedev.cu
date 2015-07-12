/*file: lebedev.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "lebedev.h"
}

#include "cuda_helper.cuh"

#define THREADS_PER_BLOCK 64

__global__ void lebedev_to_cartesian_kernel(const unsigned int lorder, double* coord)
{
  double t_a, p_a;
  double sin_t, cos_t, sin_p, cos_p;

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  if (i < lorder)
  {
    t_a = coord[i + lorder];
    p_a = coord[i + lorder * 2];

    sincos(t_a, &sin_t, &cos_t);
    sincos(p_a, &sin_p, &cos_p);

    coord[i] = sin_t * cos_p;
    coord[i + lorder] = sin_t * sin_p;
    coord[i + lorder * 2] = cos_t;
  }
}

void lebedev_to_cartesian_cuda(const unsigned int lorder, double* coord)
{
  int gridDim = ((lorder + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  lebedev_to_cartesian_kernel <<<gridDim, THREADS_PER_BLOCK>>> (lorder, coord);

  /* check if the kernel launch was successful */
  CUERR
}