/*file: atomic_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/atomic_grid.h"

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
