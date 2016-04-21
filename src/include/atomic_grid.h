/*file: grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef ATOMIC_GRID_H
#define ATOMIC_GRID_H

#include "angular.h"
#include "radial.h"

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

void atomic_grid_init(AtomicGrid *grid, AngularGrid *angular,
                      RadialGrid *radial);

/*
Calculates the integral over an atomic grid.

Args:
    segments (int): Number of array of size ``grid.size`` in the array ``f``
    f (double *): array with the value of the function ``F`` calculated in each
point of the grid.

Return:
    integral (double): value of the integral.
*/
void atomic_grid_integrate(AtomicGrid *grid, const int functions,
                           const int segments, const int size, double *f,
                           double *output);

#ifdef _CUDA

#include "cuda_helper.cuh"

#ifdef __CUDACC__

/*
Calculates the integral over an atomic grid. (CUDA kernel)

Args:
    segments (int): Number of array of size ``grid.size`` in the array ``f``
    f (double *): array with the value of the function ``F`` calculated in each
point of the grid.

Return:
    integral (double): value of the integral.
*/
__global__ void atomic_grid_integrate_kernel(const int size, const int segments,
                                             double *work, double *weights,
                                             double *integral);

/*
Calculate the multi product of a segmented array on CUDA devices.

Args:
    size (int): size of each segment.
    segments (int): Number of appended arrays on ``f``.
    f (double *): array with all the data.

Return:
    output (double *): array of size ``size`` with the multiplied arrays.
*/
__global__ void multiply_segmented_array(const int size, const int segments,
                                         double *f, double *output);
#endif
#endif

#endif