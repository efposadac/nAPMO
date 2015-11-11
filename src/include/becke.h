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

/*
CUDA functions
*/

#ifdef _CUDA

#include "cuda_helper.cuh"

struct _becke_cuda {
  int2 gridDim;   // Dim of the grid (radial, angular)
  double *radii;  // covalent radius for each center
  double *origin; // origin for each center
  double2 *xy;    // coord x and y.
  double2 *zw;    // coord z and weights .
};

typedef struct _becke_cuda BeckeGridCUDA;

/*
Copy the host grid structure into the device.
*/
void becke_init(BeckeGrid *grid, BeckeGridCUDA *grid_d);

/*
Free the memory used on the CUDA device.
*/
void becke_free(BeckeGridCUDA *grid);

#ifdef __CUDACC__

/*
Calculates the becke weights for a given ``BeckeGrid`` grid
as described in eq. 22 Becke, 1988. (CUDA Version)

References:
    Becke, A. D. A multicenter numerical integration scheme for polyatomic
molecules. J. Chem. Phys. 88, 2547 (1988).

*/
__global__ void becke_weights_kernel(BeckeGridCUDA grid, double *R_ij,
                                     double *a_ij, double *weights);

#endif

#endif
#endif
