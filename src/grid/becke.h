/*
file: becke.h
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#ifndef BECKE_H
#define BECKE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../horton/cubic_spline.h"
#include "../horton/extrapolation.h"
#include "../utils/eigen_helper.h"
#include "../utils/utils.h"
#include "atomic_grid.h"

struct BeckeGrid {

private:
  unsigned int ncenter; // Number of atoms.
  unsigned int size;    // Number of points.

  double *radii;         // covalent radius for each center
  double *origin;        // origin for each center
  double *points;        // points of the grid
  double *weights;       // weights of the grid
  double *becke_weights; // Becke weights
  double abldep;         // linear dependency tolerance por grid-based basis
  double ablmax;         // maximum l for grid-based basis

  /*
  Calculates the becke weights for a given ``BeckeGrid`` grid
  as described in eq. 22 Becke, 1988.

  References:
      Becke, A. D. A multicenter numerical integration scheme for polyatomic
  molecules. J. Chem. Phys. 88, 2547 (1988).
  */
  void compute_weights();

public:
  std::vector<AtomicGrid> atgrid;

  BeckeGrid(AtomicGrid **grids, const int n, const int l, const double ld);

  BeckeGrid(double *p, const int sz, const int nc);

  BeckeGrid(const BeckeGrid &) = default;

  BeckeGrid();

  ~BeckeGrid() {
    delete[] radii;
    delete[] origin;
    delete[] points;
    delete[] weights;
    delete[] becke_weights;
  };

  double integrate(double *f);

  double integrate(Array1D &f);

  // double *poisson_solver(double *rho, int lmax);

  int get_ncenter() { return ncenter; };

  int get_size() { return size; };

  int get_ablmax() { return ablmax; };

  double get_abldep() { return abldep; };

  double *get_radii() { return radii; };

  double *get_points() { return points; };

  double *get_weights() { return weights; };

  double *get_origin() { return origin; };

  double *get_becke_weights() { return becke_weights; };
};

#ifdef __cplusplus
extern "C" {
#endif

BeckeGrid *BeckeGrid_new(AtomicGrid **grids, int n, int l, double ld);

BeckeGrid *BeckeGrid_from_points(double *p, int sz, int nc);

void BeckeGrid_del(BeckeGrid *grid);

double BeckeGrid_integrate(BeckeGrid *grid, double *f);

int BeckeGrid_get_ncenter(BeckeGrid *grid);

int BeckeGrid_get_size(BeckeGrid *grid);

double *BeckeGrid_get_radii(BeckeGrid *grid);

double *BeckeGrid_get_points(BeckeGrid *grid);

double *BeckeGrid_get_weights(BeckeGrid *grid);

double *BeckeGrid_get_origin(BeckeGrid *grid);

double *BeckeGrid_get_becke_weights(BeckeGrid *grid);

#ifdef __cplusplus
}
#endif
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
