/*file: becke_grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef BECKE_H
#define BECKE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"

#ifdef _OMP
#include <omp.h>
#endif

struct _becke {
  int ncenter;     // Number of atoms.
  int size;        // Number of points.
  double *radii;   // covalent radius for each center
  double *origin;  // origin for each center
  double *points;  // points of the grid
  double *weights; // weights of the grid
};
typedef struct _becke BeckeGrid;

/*
Calculates the becke weights for a given ``BeckeGrid`` grid 
as described in eq. 22 Becke, 1988.

References:
    Becke, A. D. A multicenter numerical integration scheme for polyatomic
molecules. J. Chem. Phys. 88, 2547 (1988).

*/
void becke_weights(BeckeGrid *grid, double *weights);

#endif
