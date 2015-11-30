/*file: angular.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/angular.h"

void angular_cartesian(AngularGrid *grid) {
  int i, idx;
  double *x, *y, *z;

  x = (double *)malloc(grid->lorder * sizeof(double));
  y = (double *)malloc(grid->lorder * sizeof(double));
  z = (double *)malloc(grid->lorder * sizeof(double));

  ld_by_order(grid->lorder, x, y, z, grid->weights);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, idx)
#endif
  for (i = 0; i < grid->lorder; ++i) {
    idx = i * 3;
    grid->points[idx + 0] = x[i];
    grid->points[idx + 1] = y[i];
    grid->points[idx + 2] = z[i];
  }

  free(x);
  free(y);
  free(z);
}

void angular_to_spherical(AngularGrid *grid) {
  int i, idx;
  double t, p;

  // #ifdef _OMP
  // #pragma omp parallel for default(shared) private(i, idx, t, p)
  // #endif
  for (i = 0; i < grid->lorder; ++i) {
    idx = i * 3;
    xyz_to_tp(grid->points[idx + 0], grid->points[idx + 1],
              grid->points[idx + 2], &t, &p);

    grid->points[idx + 0] = 1.0;
    grid->points[idx + 1] = t;
    grid->points[idx + 2] = p;
  }
}

void angular_spherical_expansion(AngularGrid *grid, const int lmax,
                                 const int size_f, double *f, double *output) {
  int i, j, l, m, idx, lsize, lorder;
  int *lindex;
  double *t, *p, *s;

  lorder = grid->lorder;

  // Obtain Lebedev's points in spherical
  t = (double *)malloc(lorder * sizeof(double));
  p = (double *)malloc(lorder * sizeof(double));

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, idx)
#endif
  for (i = 0; i < lorder; ++i) {
    idx = i * 3;
    xyz_to_tp(grid->points[idx + 0], grid->points[idx + 1],
              grid->points[idx + 2], &t[i], &p[i]);
  }

  // Calculate index for expansion.
  lsize = (lmax + 1) * (lmax + 1);
  lindex = (int *)malloc(lsize * 2 * sizeof(int));
  idx = 0;
  for (l = 0; l <= lmax; ++l) {
    for (m = -l; m <= l; ++m) {
      lindex[idx] = l;
      lindex[lsize + idx] = m;
      idx++;
    }
  }

// Calculate expansion
#ifdef _OMP
#pragma omp parallel default(shared) private(s, i, j)
#endif
  {
    s = (double *)malloc(lorder * sizeof(double));
#ifdef _OMP
#pragma omp for
#endif
    for (i = 0; i < lsize; ++i) {
      spherical_harmonics_real(lindex[i], lindex[lsize + i], t, p, s, lorder);
      for (j = 0; j < size_f; ++j) {
        output[j * lsize + i] = angular_integrate(grid, 1, &f[j * lorder], s);
      }
    }
    free(s);
  }
  free(t);
  free(p);
  free(lindex);
}

void angular_eval_expansion(AngularGrid *grid, const int lmax, const int size_f,
                            double *decomposition, double *output) {
  int i, j, l, m, idx, lsize, lorder;
  int *lindex;
  double *t, *p, *s, aux;

  lorder = grid->lorder;

  // Obtain Lebedev's points in spherical
  t = (double *)malloc(lorder * sizeof(double));
  p = (double *)malloc(lorder * sizeof(double));

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, idx)
#endif
  for (i = 0; i < lorder; ++i) {
    idx = i * 3;
    xyz_to_tp(grid->points[idx + 0], grid->points[idx + 1],
              grid->points[idx + 2], &t[i], &p[i]);
  }

  // Calculate index for expansion.
  lsize = (lmax + 1) * (lmax + 1);
  lindex = (int *)malloc(lsize * 2 * sizeof(int));
  idx = 0;
  for (l = 0; l <= lmax; ++l) {
    for (m = -l; m <= l; ++m) {
      lindex[idx] = l;
      lindex[lsize + idx] = m;
      idx++;
    }
  }

  // Calculate spherical harmonics
  lsize = (lmax + 1) * (lmax + 1);
  s = (double *)malloc(lorder * lsize * sizeof(double));

#ifdef _OMP
#pragma omp parallel for default(shared) private(i)
#endif
  for (i = 0; i < lsize; ++i) {
    spherical_harmonics_real(lindex[i], lindex[lsize + i], t, p, &s[i * lorder], lorder);
  }

  // Evaluate expansion
#ifdef _OMP
#pragma omp parallel for default(shared) private(i, j, l, aux) collapse(2)
#endif
  for (i = 0; i < size_f; ++i) {
    for (j = 0; j < lorder; ++j) {
      aux = 0.0;
      for (l = 0; l < lsize; ++l) {
        aux += decomposition[i * lsize + l] * s[l * lorder + j];
      }
      output[i * lorder + j] = aux;
    }
  }

  free(s);
  free(t);
  free(p);
  free(lindex);
}

double angular_integrate(AngularGrid *grid, const int segments, double *f,
                         double *s) {
  int i, lorder;
  double *work;
  double output;

  lorder = grid->lorder;

  if (segments > 1) {
    work = (double *)malloc(lorder * sizeof(double));
    utils_multiply_segmented_array(lorder, segments, f, work);
  } else {
    work = &f[0];
  }

  output = 0.0;

  for (i = 0; i < lorder; ++i) {
    output += work[i] * grid->weights[i] * s[i];
  }

  if (segments > 1) {
    free(work);
  }

  return output * 4.0 * M_PI;
}
