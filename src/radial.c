/*file: radial.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/radial.h"

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
