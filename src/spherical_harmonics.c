/*file: spherical_harmonics.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "spherical_harmonics.h"

double real_spherical(int l, int m, double theta, double phi) {
  double val;

  double x = cos(theta);
  double p = gsl_sf_legendre_sphPlm(l, abs(m), x);

  if (m > 0) {
    val = sqrt(2.0) * p * cos(m * phi);
  } else if (m < 0) {
    val = sqrt(2.0) * p * sin(abs(m) * phi);
  } else {
    val = p;
  }
  return val;
};