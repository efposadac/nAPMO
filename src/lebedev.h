/*file: lebedev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef LEBEDEV_H
#define LEBEDEV_H

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*
returns the lorder number of abscissas and weights of the Lebedev quadrature
*/
void lebedev(int lorder, double * t, double* p, double* w);

#ifdef _CUDA
/*
Convert angular quadrature from spherical to cartesian coordinates (on device)
*/
void lebedev_to_cartesian_cuda(int unsigned lorder, double * coord);

#endif

#endif