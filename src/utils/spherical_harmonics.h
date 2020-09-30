/*file: spherical_harmonics.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/
#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_vector.h>

/*
Calculates the real spherical harmonic lm at theta, phi point.
Uses GNU Scientific Library GSL
*/
#ifdef __cplusplus
extern "C" {
#endif

void spherical_harmonics_real(int l, int m, double *theta, double *phi,
                              double *output, int size);

#ifdef __cplusplus
}
#endif

#endif
