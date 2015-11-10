/*file: becke_grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef BECKE_GRID_H
#define BECKE_GRID_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef _OMP
#include <omp.h>
#endif

#include "becke.h"
#include "gauss_chebyshev.h"
#include "angular.h"
#include "system.h"
  

/*
CUDA functions
*/

#ifdef _CUDA

#include "cuda_helper.cuh"

struct _grid_cuda {
  int2 gridDim; // Dim of the grid (radial, angular)
  double *radii;   // covalent radius for each center
  double *origin;  // origin for each center
  double2 *xy;  // coord x and y.
  double2 *zw;  // coord z and weights .
};

typedef struct _grid_cuda GridCuda;

/*
Copy the host grid structure into the device.
*/
void grid_init_cuda(BeckeGrid *grid, GridCuda *grid_d);

/*
Free the memory used on the CUDA device.
*/
void grid_free_cuda(GridCuda *grid);

/*
Perform the initialization (allocation and copy) of all structures needed to
perform the integration in CUDA device. Also loads additional information such
as the density matrix.
*/
double grid_integrate_cuda(System *sys, BeckeGrid *grid, double *rad,
                           int nrad);

/*
Computes the Becke weights :math:`w(r)` at point ``r`` for particle
``particleID`` as described in eq. 22 Becke, 1988. (CUDA Device Version)

References:
    Becke, A. D. A multicenter numerical integration scheme for polyatomic
molecules. J. Chem. Phys. 88, 2547 (1988).

Args:
    (double[3]): Point of the grid in which the weight will be calculated.
    particleID (int): The particle index who owns the ``r`` point.
    sys (System): system structure.

Returns:
    output (double): The value of cell_function (eq. 13, Becke, 1988) at point
``r``
*/

#ifdef __CUDACC__

__device__ double grid_weights_cuda(int n_particles, double *particle_origin,
                                    double *particle_radii, double r[3],
                                    int particleID);

/*
Functional :math:`\rho({\\bf r})`
*/
__device__ double grid_density_cuda(BasisSet basis, double *r, double *dens);

/*
Integration of the functional :math:`\rho({\\bf r})` (CUDA Version)

Kernel to perform the integration based on multicenter molecular integration
from Becke.
Please note that the structures are not passed by reference but value. All
structure members have to be allocated/copied properly before to use this
kernel.

TODO:
    Pass function pointer to enable the use of different kind of functionals.
*/
__global__ void
grid_integrate_kernel(const System sys, const int2 gridDim,
                      double *__restrict__ radii, double *__restrict__ origin,
                      double2 *__restrict__ xy, double2 *__restrict__ zw,
                      double *__restrict__ dens, double *__restrict__ rad,
                      int nrad, double *__restrict__ integral);

/*
Iterated cutoff profile. eq. 21, Becke 1988. (CUDA Device version)

Note:
    Avoid the use of __host__ __device__ in order to allow the compilation with
    other compilers in the case of non CUDA compilation.

*/
__device__ __forceinline__ double grid_soft_mu_cuda(const double &mu) {
  return 0.5 * mu * (3.0 - mu * mu);
}

#endif

#endif

#endif
