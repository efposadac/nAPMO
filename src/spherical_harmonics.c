/*file: spherical_harmonics.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "spherical_harmonics.h"

void real_spherical(int l, int m, double *theta, double *phi, double* output, int size) {
  
  int i;
  double x;
  double p;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, x, p)
#endif
  for (i = 0; i < size; ++i)
  {
    x = cos(theta[i]);
    p = gsl_sf_legendre_sphPlm(l, abs(m), x);  
    if (m > 0) {
      output[i] = sqrt(2.0) * p * cos(m * phi[i]);
    } else if (m < 0) {
      output[i] = sqrt(2.0) * p * sin(abs(m) * phi[i]);
    } else {
      output[i] = p;
    }
  }
};