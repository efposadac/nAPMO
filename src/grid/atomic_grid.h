/*file: grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef ATOMIC_GRID_H
#define ATOMIC_GRID_H

#include <iostream>

#include "angular.h"
#include "radial.h"

#include "../utils/omp_helper.h"

struct AtomicGrid {

private:
  unsigned int size; // Number of points.
  double radii;      // covalent radius for center
  double *origin;    // origin for each center
  double *points;    // points of the grid
  double *weights;   // weights of the grid

public:
  AtomicGrid() = default;

  AtomicGrid(const AtomicGrid &) = default;

  AtomicGrid(AngularGrid *angular, RadialGrid *radial, double *R);

  ~AtomicGrid() {
    delete[] origin;
    delete[] points;
    delete[] weights;
  };
  /*
  Calculates the integral over an atomic grid.

  Args:
      segments (int): Number of array of size ``grid.size`` in the array ``f``
      f (double *): array with the value of the function ``F`` calculated in
                    each point of the grid.

  Return:
      integral (double): value of the integral.
  */
  double *integrate(const unsigned int nfunc, const unsigned int segments,
                    const unsigned int size, double *f);

  unsigned int get_size() { return size; };
  double get_radii() { return radii; };
  double *get_origin() { return origin; };
  double *get_points() { return points; };
  double *get_weights() { return weights; };
};

#ifdef __cplusplus
extern "C" {
#endif

AtomicGrid *AtomicGrid_new(AngularGrid *angular, RadialGrid *radial, double *R);

void AtomicGrid_del(AtomicGrid *grid);

double *AtomicGrid_integrate(AtomicGrid *grid, int nfunc, int segments,
                             int size, double *f);

int AtomicGrid_get_size(AtomicGrid *grid);

double AtomicGrid_get_radii(AtomicGrid *grid);

double *AtomicGrid_get_origin(AtomicGrid *grid);

double *AtomicGrid_get_points(AtomicGrid *grid);

double *AtomicGrid_get_weights(AtomicGrid *grid);

#ifdef __cplusplus
}
#endif

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