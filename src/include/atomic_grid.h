/*file: grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef ATOMIC_GRID_H
#define ATOMIC_GRID_H

#include "angular.h"

#ifdef _OMP
#include <omp.h>
#endif

struct _atomic_grid {
  int size;        // Number of points.
  double radii;    // covalent radius for each center
  double *origin;  // origin for each center
  double *points;  // points of the grid
  double *weights; // weights of the grid
};
typedef struct _atomic_grid AtomicGrid;

double atomic_grid_integrate(AtomicGrid *grid, const int segments, double *f);

#ifdef _CUDA

#include "cuda_helper.cuh"

#ifdef __CUDACC__
__global__ void atomic_grid_integrate_kernel(const int size, const int segments,
                                             double *work, double *weights,
                                             double *integral);

__global__ void multiply_segmented_array(const int size, const int segments,
                                         double *f, double *output);
#endif
#endif

#endif