/*file: becke_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

//#define _GNU_SOURCE
#include "include/becke_grid.h"

double grid_weights(BeckeGrid *grid, double r[3], int particleID) {
  double x_i[3], x_j, r_i, r_j, R_ij, mu_ij, rm_i, rm_j, chi, u_ij, a_ij, nu_ij;
  double sum = 0, output, aux1, aux2, aux3 = 0;
  int ncenter = grid->ncenter;
  int i, j, k;

  for (i = 0; i < ncenter; ++i) {
    aux1 = 1.0;
    x_i[0] = grid->origin[i * 3 + 0];
    x_i[1] = grid->origin[i * 3 + 1];
    x_i[2] = grid->origin[i * 3 + 2];
    rm_i = grid->radii[i];

    for (j = 0; j < ncenter; ++j) {
      if (i != j) {
        // Internuclear distance (R_ij eq. 11)
        r_i = r_j = R_ij = 0.0;

        for (k = 0; k < 3; ++k) {
          x_j = grid->origin[j * 3 + k];
          r_i += (r[k] - x_i[k]) * (r[k] - x_i[k]);
          r_j += (r[k] - x_j) * (r[k] - x_j);
          R_ij += (x_i[k] - x_j) * (x_i[k] - x_j);
        }

        r_i = sqrt(r_i);
        r_j = sqrt(r_j);
        R_ij = sqrt(R_ij);

        // \mu_ij eq. 11
        mu_ij = (r_i - r_j) / R_ij;

        // Atomic size adjustment. see appendix, Becke, 1988.
        rm_j = grid->radii[j];

        // eq. A4
        chi = rm_i / rm_j;
        // eq. A6
        u_ij = (chi - 1.0) / (chi + 1.0);
        // eq. A5
        a_ij = u_ij / ((u_ij * u_ij) - 1.0);
        // eq. A3
        if (fabs(a_ij) > 0.50) {
          a_ij *= 0.50 / fabs(a_ij);
        }
        // eq. A2
        nu_ij = mu_ij + a_ij * (1.0 - mu_ij * mu_ij);

        aux2 = grid_soft_mu(grid_soft_mu(grid_soft_mu(nu_ij)));
        aux1 *= 0.5 * (1.0 - aux2);
      }
    }
    sum += aux1;
    if (i == particleID)
      aux3 = aux1;
  }
  // eq. 22
  output = aux3 / sum;
  return output;
}

double grid_integrate(BeckeGrid *grid, System *sys, double *rad, int nrad) {

#ifdef _CUDA
/*
TODO: Change this call for a "try:" style. The idea is that if there is not
device available for CUDA the code will be executed on the host.
*/
return grid_integrate_cuda(sys, grid, rad, nrad);

#endif
  double integral, p;
  int i, j, idx, idxr, size, size2;

  int ncenter = grid->ncenter;

  // Fetch density file.
  size = sys->basis.n_cont * sys->basis.n_cont;
  double *dens = (double *)malloc(size * sizeof(double));

  FILE *file;
  file = fopen("data.dens", "r");

  for (i = 0; i < size; ++i) {
    if (!fscanf(file, "%lf", &dens[i])) {
      exit(1);
    };
  }

  integral = 0.0;
  idxr = 0;
  size = grid->size / ncenter;
  size2 = size / nrad;
#ifdef _OMP
#pragma omp parallel for default(shared) private(i, j, idx, idxr,              \
                                                 p) reduction(+ : integral)
#endif
  for (i = 0; i < ncenter; ++i) {
    for (j = 0; j < size; ++j) {

      // Calculate Becke weights
      idx = (i * size + j) * 3;
      p = grid_weights(grid, &grid->points[idx], i);

      // Calculate Integral
      idxr = i * nrad + floor(j / size2);
      integral +=
          (rad[idxr] * rad[idxr] * p * grid->weights[j] * grid->radii[i] *
           grid_density(sys, &grid->points[idx], dens));
    }
  }

  free(dens);
  return integral * 4.0 * M_PI;
}

/*
Temporal functional (test)
*/
double grid_density(System *sys, double *r, double *dens) {
  BasisSet basis = sys->basis;
  int i, j, aux, counter = 0;
  int n_cont = basis.n_cont;
  double temp;
  double function_value, output = 0.0;
  double RP2, factor;

  double *basis_val = (double *)calloc(n_cont * 2, sizeof(double));

  for (i = 0; i < n_cont; ++i) {
    factor = 1.0, RP2 = 0.0;
    for (j = 0; j < 3; ++j) {
      aux = i * 3;
      temp = r[j] - basis.origin[aux + j];
      RP2 += (temp * temp);
      factor *= pow(temp, basis.basis_l[aux + j]);
    }

    function_value = 0.0;
    for (j = 0; j < basis.n_prim_cont[i]; ++j) {
      function_value +=
          basis.coefficient[counter] * exp(-basis.exponent[counter] * RP2);
      counter += 1;
    }

    function_value *= factor * basis.normalization[i];

    for (j = 0; j < n_cont; ++j) {
      basis_val[n_cont + j] += function_value * dens[i * n_cont + j];
    }

    basis_val[i] = function_value;
  }

  for (i = 0; i < n_cont; ++i) {
    output += basis_val[i] * basis_val[n_cont + i];
  }

  free(basis_val);
  return output;
}
