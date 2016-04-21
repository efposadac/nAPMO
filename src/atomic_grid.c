/*file: atomic_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "include/atomic_grid.h"

void atomic_grid_init(AtomicGrid *grid, AngularGrid *angular,
                      RadialGrid *radial) {

  int i, j, idx, idy, aux;
  double ap[3], aw, rp, rw;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, j, aux, idx, idy, ap, aw,  \
                                                 rp, rw)
#endif
  for (i = 0; i < radial->size; ++i) {
    rp = radial->points[i];
    rw = radial->weights[i];

    for (j = 0; j < angular->lorder; ++j) {
      idx = j * 3;
      ap[0] = angular->points[idx + 0];
      ap[1] = angular->points[idx + 1];
      ap[2] = angular->points[idx + 2];

      aw = angular->weights[j];

      aux = i * angular->lorder + j;
      idy = aux * 3;

      grid->points[idy + 0] = ap[0] * rp + grid->origin[0];
      grid->points[idy + 1] = ap[1] * rp + grid->origin[1];
      grid->points[idy + 2] = ap[2] * rp + grid->origin[2];

      grid->weights[aux] = aw * rw;
    }
  }
}

#ifndef _CUDA
void atomic_grid_integrate(AtomicGrid *grid, const int functions,
                           const int segments, const int size, double *f,
                           double *output) {

  int i, j;
  int total_size = size * segments;
  double *work;

  // Convert in one array of size ``total_size``
  if (functions > 1) {
    work = (double *)malloc(total_size * sizeof(double));
    utils_multiply_segmented_array(total_size, functions, f, work);
  } else {
    work = &f[0];
  }

  //OMP?
  
  // Perform reduction
  for (i = 0; i < segments; ++i) {
    for (j = 0; j < size; ++j) {
      output[i] += work[i*size+j];
    }
  }

  if (functions > 1) {
    free(work);
  }
}

#endif