/*file: lebedev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "lebedev.h"

void lebedev_spherical(const int lorder, double *t, double *p, double *w) {
  int i;
  double *x, *y, *z;

  x = (double *)malloc(lorder * sizeof(double));
  y = (double *)malloc(lorder * sizeof(double));
  z = (double *)malloc(lorder * sizeof(double));

  ld_by_order(lorder, x, y, z, w);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < lorder; ++i) {
    xyz_to_tp(x[i], y[i], z[i], &t[i], &p[i]);
  }

  free(x);
  free(y);
  free(z);
}

void lebedev_cartesian(const int lorder, double *points, double *w) {
  int i, idx;
  double *x, *y, *z;

  x = (double *)malloc(lorder * sizeof(double));
  y = (double *)malloc(lorder * sizeof(double));
  z = (double *)malloc(lorder * sizeof(double));

  ld_by_order(lorder, x, y, z, w);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, idx)
#endif
  for (i = 0; i < lorder; ++i) {
    idx = i * 3;
    points[idx + 0] = x[i];
    points[idx + 1] = y[i];
    points[idx + 2] = z[i];
  }

  free(x);
  free(y);
  free(z);
}

void lebedev_spherical_expansion(const int lorder, const int lmax,
                                 const int size, double *f, double *output) {
  int i, l, m, idx, lsize;
  double *t, *p, *w, *s;

  // Obtain Lebedev's points
  t = (double *)malloc(lorder * sizeof(double));
  p = (double *)malloc(lorder * sizeof(double));
  w = (double *)malloc(lorder * sizeof(double));

  lebedev_spherical(lorder, t, p, w);

  // Calculate expansion
  s = (double *)malloc(lorder * 2 * sizeof(double));

  idx = 0;
  lsize = (lmax + 1) * (lmax + 1);
  for (l = 0; l <= lmax; ++l) {
    for (m = -l; m <= l; ++m) {
      real_spherical(l, m, t, p, s, lorder);
      for (i = 0; i < size; ++i) {
        memcpy(&s[lorder], &f[i * lorder], lorder * sizeof(double));
        output[i * lsize + idx] = lebedev_integrate(lorder, 2, s, w);
      }
      idx++;
    }
  }

  free(t);
  free(p);
  free(w);
  free(s);
}

double lebedev_integrate(const int lorder, const int nfunc, double *f,
                         double *w) {
  int i;
  double *work;
  double output;

  work = (double *)malloc(lorder * sizeof(double));
  utils_multiply_segmented_array(lorder, nfunc, f, work);

  output = 0.0;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i) reduction(+ : output)
#endif
  for (i = 0; i < lorder; ++i) {
    output += work[i] * w[i];
  }

  free(work);
  return output * 4.0 * M_PI;
}
