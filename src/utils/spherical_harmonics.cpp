/*file: spherical_harmonics.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

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

std::complex<double> spherical_harmonics_complex(int l, int m, double theta,
                                                 double phi) {

  // printf("spherical_harmonics: l: %d m: %d\n", l, m);
  
  std::complex<double> result(cos(abs(m) * phi), sin(abs(m) * phi));

  double x = cos(theta);
  double p = gsl_sf_legendre_sphPlm(l, abs(m), x);

  result *= p;

  if (m < 0.0) {
    double sign = (m % 2) ? -1.0 : 1.0;
    result = sign * std::conj(result);
  }

  return result;
}
