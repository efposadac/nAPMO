/*file: poisson_helper.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef POISSON_HELPER_H
#define POISSON_HELPER_H

#include "radial.h"

void finite_difference_matrix(RadialGrid *rgrid, double *A, int l);

#endif