/*file: becke.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

extern "C" {
#include "include/becke.h"
}

#define THREADS_PER_BLOCK 64

/*
Implementation of functions, for documentation see the header file.
*/

void becke_init(BeckeGrid *grid, BeckeGridCUDA *grid_d) {
  grid_d->gridDim = make_int2(grid->size, grid->ncenter);

  /* Allocate space for grid on device*/
  int ncenter = grid_d->gridDim.y * sizeof(double);
  cudaMalloc((void **)&grid_d->radii, ncenter);
  cudaMalloc((void **)&grid_d->origin, ncenter * 3);

  int size = grid_d->gridDim.x * sizeof(double2);
  cudaMalloc((void **)&grid_d->xy, size);
  cudaMalloc((void **)&grid_d->zw, size);

  /* Copying atomic information*/
  {
    cudaMemcpy(grid_d->radii, grid->radii, ncenter, cudaMemcpyHostToDevice);
    cudaMemcpy(grid_d->origin, grid->origin, ncenter * 3,
               cudaMemcpyHostToDevice);
  }

  /* Copying grid points*/
  {
    double2 *buffer_xy = (double2 *)malloc(size);
    double2 *buffer_zw = (double2 *)malloc(size);

    for (int i = 0; i < grid_d->gridDim.x; ++i) {
      int idx = i * 3;
      buffer_xy[i] = make_double2(grid->points[idx], grid->points[idx + 1]);
      buffer_zw[i] = make_double2(grid->points[idx + 2], grid->weights[i]);
    }

    cudaMemcpy(grid_d->xy, buffer_xy, size, cudaMemcpyHostToDevice);
    cudaMemcpy(grid_d->zw, buffer_zw, size, cudaMemcpyHostToDevice);

    free(buffer_xy);
    free(buffer_zw);
  }
}

void becke_free(BeckeGridCUDA *grid) {
  cudaFree(grid->radii);
  cudaFree(grid->origin);
  cudaFree(grid->xy);
  cudaFree(grid->zw);
}

void becke_weights(BeckeGrid *grid, double *weights) {
  int iatom, jatom, npoint, offset, size;
  double *R_ij, *a_ij, aux;

  size = (grid->ncenter * (grid->ncenter + 1)) / 2;

  // Calculate internuclear distance
  R_ij = (double *)malloc(size * sizeof(double));

  offset = 0;
  for (iatom = 0; iatom < grid->ncenter; iatom++) {
    for (jatom = 0; jatom <= iatom; jatom++) {
      R_ij[offset] =
          utils_distance(&grid->origin[3 * iatom], &grid->origin[3 * jatom]);
      offset += 1;
    }
  }

  // Calculate atomic size parameter
  a_ij = (double *)malloc(size * sizeof(double));
  offset = 0;
  for (iatom = 0; iatom < grid->ncenter; iatom++) {
    for (jatom = 0; jatom <= iatom; jatom++) {
      aux = grid->radii[iatom] / grid->radii[jatom];
      // eq. A6
      aux = (aux - 1.0) / (aux + 1.0);
      // eq. A5
      aux = aux / ((aux * aux) - 1.0);
      // eq. A3
      if (fabs(aux) > 0.50) {
        aux *= 0.50 / fabs(aux);
      }

      a_ij[offset] = aux;
      offset += 1;
    }
  }

  // Prepare Cuda data
  double *R_ij_d, *a_ij_d, *weights_d;

  cudaMalloc((void **)&R_ij_d, size * sizeof(double));
  cudaMalloc((void **)&a_ij_d, size * sizeof(double));
  cudaMalloc((void **)&weights_d, grid->size * sizeof(double));

  cudaMemcpy(R_ij_d, R_ij, size * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(a_ij_d, a_ij, size * sizeof(double), cudaMemcpyHostToDevice);

  BeckeGridCUDA grid_d;
  becke_init(grid, &grid_d);

  // Computation of Becke weights
  npoint = grid_d.gridDim.x / grid_d.gridDim.y;
  dim3 dimGrid(((npoint + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1, 1);

  becke_weights_kernel<<<dimGrid, THREADS_PER_BLOCK>>>(grid_d, R_ij_d, a_ij_d,
                                                       weights_d);
  CUERR

  /*
  Bring result back
  */
  cudaMemcpy(weights, weights_d, grid->size * sizeof(double),
             cudaMemcpyDeviceToHost);

  free(R_ij);
  free(a_ij);
  cudaFree(R_ij_d);
  cudaFree(a_ij_d);
  cudaFree(weights_d);
  becke_free(&grid_d);
}

__global__ void becke_weights_kernel(BeckeGridCUDA grid, double *R_ij,
                                     double *a_ij, double *weights) {

  const unsigned int point = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  const unsigned int npoint = grid.gridDim.x / grid.gridDim.y;

  int offset, order, idx;
  double r[3], r_i, r_j, mu_ij;
  double sum, aux, p, s;

  if (point < npoint) {
    order = 3;
    for (int atom = 0; atom < grid.gridDim.y; ++atom) {

      idx = atom * npoint + point;
      const double2 aux_xy = grid.xy[idx];
      const double2 aux_zw = grid.zw[idx];

      r[0] = aux_xy.x;
      r[1] = aux_xy.y;
      r[2] = aux_zw.x;

      sum = 0.0;
      aux = 0.0;
      for (int iatom = 0; iatom < grid.gridDim.y; ++iatom) {

        p = 1.0;
        for (int jatom = 0; jatom < grid.gridDim.y; ++jatom) {

          if (iatom == jatom)
            continue;

          // offset
          if (iatom < jatom) {
            offset = (jatom * (jatom + 1)) / 2 + iatom;
          } else {
            offset = (iatom * (iatom + 1)) / 2 + jatom;
          }

          // Internuclear distance (R_ij eq. 11)
          r_i = distance(r, &grid.origin[iatom * 3]);
          r_j = distance(r, &grid.origin[jatom * 3]);

          // \mu_ij eq. 11
          mu_ij = (r_i - r_j) / R_ij[offset];

          // eq. A2
          s = mu_ij + a_ij[offset] * (1.0 - mu_ij * mu_ij);

          // eq. 19 and 20
          for (int k = 1; k <= order; k++) {
            s = 0.5 * s * (3 - s * s);
          }

          s = 0.5 * (1.0 - s);
          p *= s;
        }

        sum += p;
        if (iatom == atom)
          aux = p;
      }
      // eq. 22
      weights[idx] = aux / sum;
    }
  }
}
