/*file: grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "grid.h"

void grid_evaluate_atomic_expansion(int lmax, int lorder, int size,
                                        double *decomposition, double *points,
                                        double *output) {
  int i, j, l, m, idx, lsize;
  double *t, *p, *w, *s;

  // Obtain Lebedev's points
  t = (double *)malloc(lorder * sizeof(double));
  p = (double *)malloc(lorder * sizeof(double));
  w = (double *)malloc(lorder * sizeof(double));

  lebedev_spherical(lorder, t, p, w);

  // Calculate spherical harmonics
  lsize = (lmax + 1) * (lmax + 1);
  s = (double *)malloc(lorder * lsize * sizeof(double));

  idx = 0;
  for (l = 0; l <= lmax; ++l) {
    for (m = -l; m <= l; ++m) {
      real_spherical(l, m, t, p, &s[idx * lorder], lorder);
      idx++;
    }
  }

  // Evaluate expansion
  idx = 0;
  for (i = 0; i < size; ++i) {
    for (j = 0; j < lorder; ++j) {
      output[idx] = 0.0;
      for (l = 0; l < lsize; ++l) {
        output[idx] += decomposition[i * lsize + l] * s[l * lorder + j];
      }
      idx++;
    }
  }

  free(t);
  free(p);
  free(w);
  free(s);
}

double grid_atomic_integrate(const int size, const int segments, double rm, double *f,
                         double *w){
  int i;
  double *work;
  double output;

  work = (double *)malloc(size * sizeof(double));
  multiply_segmented_array(size, segments, f, work);

  output = 0.0;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i) reduction(+ : output)
#endif
  for (i = 0; i < size; ++i) {
    output += work[i] * w[i];
  }

  free(work);
  return output * 4.0 * M_PI * rm;   
}