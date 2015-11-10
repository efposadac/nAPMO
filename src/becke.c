/*file: becke_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/becke.h"

#ifndef _CUDA

void becke_weights(BeckeGrid *grid, double *weights) {
  int atom, iatom, jatom, point, npoint, offset, size, order, k, idx;
  double r_i, r_j, mu_ij;
  double sum = 0.0, aux, p = 1.0, s;
  double *R_ij;
  double *a_ij;

  order = 3;
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

  // Computation of Becke weights
  npoint = grid->size / grid->ncenter;

  idx = 0;
  for (atom = 0; atom < grid->ncenter; ++atom) {
#ifdef _OMP
#pragma omp parallel for default(shared)                                       \
    firstprivate(atom, npoint, idx) private(point, iatom, jatom, offset, r_i,  \
                                            r_j, mu_ij, s, k, aux, sum, p)
#endif
    for (point = 0; point < npoint; ++point) {
      sum = 0.0;
      aux = 0.0;
      for (iatom = 0; iatom < grid->ncenter; ++iatom) {
        p = 1.0;
        for (jatom = 0; jatom < grid->ncenter; ++jatom) {
          if (iatom == jatom)
            continue;

          // offset
          if (iatom < jatom) {
            offset = (jatom * (jatom + 1)) / 2 + iatom;
          } else {
            offset = (iatom * (iatom + 1)) / 2 + jatom;
          }

          // Internuclear distance (R_ij eq. 11)
          idx = (atom * npoint + point) * 3;
          r_i = utils_distance(&grid->points[idx], &grid->origin[iatom * 3]);
          r_j = utils_distance(&grid->points[idx], &grid->origin[jatom * 3]);

          // \mu_ij eq. 11
          mu_ij = (r_i - r_j) / R_ij[offset];

          // eq. A2
          s = mu_ij + a_ij[offset] * (1.0 - mu_ij * mu_ij);

          // eq. 19 and 20
          for (k = 1; k <= order; k++) {
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
      weights[atom * npoint + point] = aux / sum;
    }
  }

  free(R_ij);
  free(a_ij);
}

#endif