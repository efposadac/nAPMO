/*file: poisson_helper.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/poisson_helper.h"

void finite_difference_matrix(RadialGrid *rgrid, double *data, int *row,
                              int *col, int l) {

  int i, j, idx;
  int npoints = rgrid->size;
  double dz, d2z, point;

  double h = rgrid->z[0];
  double h2_inv = 1.0 / (h * h);
  double h_inv = 1.0 / h;

  double dcoeff[3] = {-0.5, 0.0, 0.5};
  double d2coeff[3] = {1.0, -2.0, 1.0};
  double aux[3] = {0.0, 0.0, 0.0};

  data[0] = 1.0;
  row[0] = 0;
  col[0] = 0;

  aux[1] = l * (l + 1.0);

// #ifdef _OMP
// #pragma omp parallel for default(shared)                                       \
//     firstprivate(aux) private(i, j, idx, point, dz, d2z)
// #endif
  for (i = 0; i < npoints; ++i) {
    point = rgrid->points[i];

    // Build A
    dz = rgrid->dz[i] * rgrid->dz[i];
    d2z = rgrid->d2z[i];

    aux[1] /= (point * point);
    dz *= h2_inv;
    d2z *= h_inv;

    idx = i * 3 + 1;
    for (j = 0; j < 3; ++j) {
      row[idx + j] = i + 1;
      col[idx + j] = i + j;
      data[idx + j] = (d2coeff[j] * dz) + (dcoeff[j] * d2z) - aux[j];
    }
  }

  idx = npoints * 3 + 1;
  data[idx] = 1.0;
  row[idx] = npoints + 1;
  col[idx] = npoints + 1;
}
