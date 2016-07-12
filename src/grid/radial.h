/*file: radial.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef RADIAL_H
#define RADIAL_H

#include <math.h>
#include <stdlib.h>
#include "../utils/utils.h"
#include "gauss_chebyshev.h"

struct _radial {
  int size;        // number of grid points.
  double radii;    // covalent radius of atom
  double *points;  // points of the grid
  double *weights; // weights of the grid
  double *z;       // Uniform spaced grid
  double *dz;      // first derivative of z
  double *d2z;     // second derivative of z
};
typedef struct _radial RadialGrid;

/*
Returns the radial points mapped uniformly in the interval [0,1], see
Becke's paper.
*/
void radial_get_z(RadialGrid *grid);

/*
Returns the first derivative of the uniform z grid.
*/
void radial_deriv_z(RadialGrid *grid);

/*
Returns the second derivative of the uniform z grid.
*/
void radial_deriv2_z(RadialGrid *grid);

double radial_integrate(RadialGrid *grid, const int segments, double *f);

#endif