/*file: spherical_harmonics.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "spherical_harmonics.h"

void spherical_harmonics_real(int l, int m, double *theta, double *phi,
                              double *output, int size) {

  // printf("spherical_harmonics: l: %d m: %d\n", l, m);

  for (int i = 0; i < size; ++i) {
    double x = cos(theta[i]);
    double p = gsl_sf_legendre_sphPlm(l, abs(m), x);
    if (m > 0) {
      output[i] = sqrt(2.0) * p * cos(m * phi[i]);
    } else if (m < 0) {
      output[i] = sqrt(2.0) * p * sin(abs(m) * phi[i]);
    } else {
      output[i] = p;
    }
  }
};