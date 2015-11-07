/*file: grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef ATOMIC_GRID_H
#define ATOMIC_GRID_H

#include "angular.h"

struct _atomic_grid {
  int size;        // Number of points.
  double radii;    // covalent radius for each center
  double *origin;  // origin for each center
  double *points;  // points of the grid
  double *weights; // weights of the grid
};
typedef struct _atomic_grid AtomicGrid;

double atomic_grid_integrate(AtomicGrid *grid, const int segments, double *f);
#endif