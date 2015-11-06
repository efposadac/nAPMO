/*file: poisson_helper.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "poisson_helper.h"

void finite_difference_matrix(RadialGrid *rgrid, double *A, int l) {

  int i, j;
  int npoints = rgrid->size;
  int npoints2 = npoints + 2;
  double dz, d2z, aux, point;

  double h = rgrid->z[0];
  double h2_inv = 1.0 / (h * h);
  double h_inv = 1.0 / h;

  double dcoeff[3] = {-0.5, 0.0, 0.5};
  double d2coeff[3] = {1.0, -2.0, 1.0};

  A[0] = 1.0;
  A[((npoints + 2) * (npoints + 2)) - 1] = 1.0;

  aux = l * (l + 1.0);

  for (i = 0; i < npoints; ++i) {
    point = rgrid->points[i];

    // Build A
    dz = rgrid->dz[i] * rgrid->dz[i];
    d2z = rgrid->d2z[i];

    aux /= (point * point);
    dz *= h2_inv;
    d2z *= h_inv;

    for (j = 0; j < 3; ++j) {
      A[(i + 1) * npoints2 + (i + j)] = (d2coeff[j] * dz) + (dcoeff[j] * d2z);
      A[(i + 1) * npoints2 + (i + 1)] += -aux;
    }
  }
}