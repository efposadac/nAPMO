/*file: radial.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "radial.h"

#include "stdio.h"

void radial_init(RadialGrid *grid) {
  gaussChebyshev(grid->size, grid->radii, grid->points, grid->weights);
}

void radial_get_z(RadialGrid *grid) {
  int i;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < grid->size; ++i) {
    grid->z[i] = acos((grid->points[i] - grid->radii) /
                      (grid->points[i] + grid->radii)) *
                 M_1_PI;
  }
}

void radial_deriv_z(RadialGrid *grid) {
  int i;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < grid->size; ++i) {
    grid->dz[i] = -sqrt(grid->radii * grid->points[i] /
                        pow(grid->radii + grid->points[i], 2)) /
                  (M_PI * grid->points[i]);
  }
}

void radial_deriv2_z(RadialGrid *grid) {
  int i;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < grid->size; ++i) {
    grid->d2z[i] = grid->radii * grid->radii *
                   (grid->radii + (3.0 * grid->points[i])) /
                   (2.0 * M_PI * pow(grid->radii * grid->points[i] /
                                         pow(grid->radii + grid->points[i], 2),
                                     1.5) *
                    pow(grid->radii + grid->points[i], 5));
  }
}

double radial_integrate(RadialGrid *grid, const int segments, double *f) {
  int i, size;
  double *work;
  double output;

  size = grid->size;

  if (segments > 1) {
    work = (double *)malloc(size * sizeof(double));
    utils_multiply_segmented_array(size, segments, f, work);
  } else {
    work = &f[0];
  }

  output = 0.0;

  for (i = 0; i < size; ++i) {
    output += work[i] * grid->weights[i];
  }

  if (segments > 1) {
    free(work);
  }

  return output;
}