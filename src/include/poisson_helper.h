/*file: poisson_helper.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef POISSON_HELPER_H
#define POISSON_HELPER_H

#include "radial.h"

/*
Build the A matrix to solve a second order PDE, eq. 21 Becke's paper.
*/
void finite_difference_matrix(RadialGrid *rgrid, double *data, int *row, int *col, int l);

#endif