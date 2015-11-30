/*file: atomic_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

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

      grid->weights[aux] = aw * rw * rp * rp;
    }
  }
}

#ifndef _CUDA
double atomic_grid_integrate(AtomicGrid *grid, const int segments, double *f) {
  int i;
  double *work;
  double output;

  int size = grid->size;

  if (segments > 1) {
    work = (double *)malloc(size * sizeof(double));
    utils_multiply_segmented_array(size, segments, f, work);
  } else {
    work = &f[0];
  }

  output = 0.0;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i) reduction(+ : output)
#endif
  for (i = 0; i < size; ++i) {
    output += work[i] * grid->weights[i];
  }

  free(work);

  return output * 4.0 * M_PI * grid->radii;
}

#endif