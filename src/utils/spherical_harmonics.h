/*file: spherical_harmonics.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/
#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_legendre.h>
#include <stdlib.h>
#include <math.h>

/*
Calculates the real spherical harmonic lm at theta, phi point.
Uses GNU Scientific Library GSL
*/
void spherical_harmonics_real(int l, int m, double *theta, double *phi, double *output,
                    int size);

#endif