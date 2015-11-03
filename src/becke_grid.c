/*file: becke_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

//#define _GNU_SOURCE
#include "becke_grid.h"

void grid_init(Grid *grid) {
  grid->radial_abscissas = (double *)malloc(grid->n_radial * sizeof(double));
  grid->radial_weights = (double *)malloc(grid->n_radial * sizeof(double));

  grid->angular_theta = (double *)malloc(grid->n_angular * sizeof(double));
  grid->angular_phi = (double *)malloc(grid->n_angular * sizeof(double));
  grid->angular_weights = (double *)malloc(grid->n_angular * sizeof(double));

  lebedev_spherical(grid->n_angular, grid->angular_theta, grid->angular_phi,
                    grid->angular_weights);
  gaussChebyshev(grid->n_radial, 1.0, grid->radial_abscissas,
                 grid->radial_weights);
}

void grid_free(Grid *grid) {
  free(grid->radial_abscissas);
  free(grid->radial_weights);
  free(grid->angular_theta);
  free(grid->angular_phi);
  free(grid->angular_weights);
}

double grid_weights(System *sys, double r[3], int particleID) {
  double x_i[3], x_j, r_i, r_j, R_ij, mu_ij, rm_i, rm_j, chi, u_ij, a_ij, nu_ij;
  double sum = 0, output, aux1, aux2, aux3 = 0;
  int n_particles = sys->n_particles;
  int i, j, k;

  for (i = 0; i < n_particles; ++i) {
    aux1 = 1.0;
    x_i[0] = sys->particle_origin[i * 3 + 0];
    x_i[1] = sys->particle_origin[i * 3 + 1];
    x_i[2] = sys->particle_origin[i * 3 + 2];
    rm_i = sys->particle_radii[i];

    for (j = 0; j < n_particles; ++j) {
      if (i != j) {
        // Internuclear distance (R_ij eq. 11)
        r_i = r_j = R_ij = 0.0;

        for (k = 0; k < 3; ++k) {
          x_j = sys->particle_origin[j * 3 + k];
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
        rm_j = sys->particle_radii[j];

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

double grid_integrate(System *sys, Grid *grid) {
#ifdef _CUDA
  /*
  TODO: Change this call for a "try:" style. The idea is that if there is not
  device available for CUDA the code will be executed on the host.
  */
  return grid_integrate_cuda(sys, grid);
#endif
  double integral, rm, rad, r[3];
  double q_r, w_r, t_a, p_a, w_a, p;
  double sin_t, cos_t, sin_p, cos_p;

  int n_radial = grid->n_radial;
  int n_angular = grid->n_angular;
  int n_particles = sys->n_particles;
  int i, j, k, aux5;

  // Fetch density file.
  int size = sys->basis.n_cont * sys->basis.n_cont;
  double *dens = (double *)malloc(size * sizeof(double));

  FILE *file;
  file = fopen("data.dens", "r");

  for (i = 0; i < size; ++i) {
    if (!fscanf(file, "%lf", &dens[i])) {
      exit(1);
    };
  }

  integral = 0.0;

#ifdef _OMP
/*
Note: I've implemented the unroll option for inner loops without any impact on
the performance.
*/
#pragma omp parallel for default(shared) private(i, j, k, q_r, w_r, t_a, p_a,  \
                                                 w_a, sin_t, cos_t, sin_p,     \
                                                 cos_p, rm, rad, aux5, r,      \
                                                 p) reduction(+ : integral)
#endif
  for (i = 0; i < n_radial; ++i) {
    q_r = grid->radial_abscissas[i];
    w_r = grid->radial_weights[i];

    for (j = 0; j < n_angular; ++j) {
      t_a = grid->angular_theta[j];
      p_a = grid->angular_phi[j];
      w_a = grid->angular_weights[j];
      sincos(t_a, &sin_t, &cos_t);
      sincos(p_a, &sin_p, &cos_p);

      for (k = 0; k < n_particles; ++k) {
        rm = sys->particle_radii[k];
        if (sys->particle_number[k] != 1) {
          rm *= 0.5;
        }

        rad = rm * q_r;
        aux5 = k * 3;

        r[0] = (rad * sin_t * cos_p) + sys->particle_origin[aux5 + 0];
        r[1] = (rad * sin_t * sin_p) + sys->particle_origin[aux5 + 1];
        r[2] = (rad * cos_t) + sys->particle_origin[aux5 + 2];

        // Calculate Becke weights
        p = grid_weights(sys, r, k);

        // Calculate Integral
        integral +=
            (rad * rad * p * w_r * w_a * rm * grid_density(sys, r, dens));
      }
    }
  }

  free(dens);
  return integral * 4.0 * M_PI;
}

/*
Temporal functional (test)
*/
double grid_density(System *sys, double *r, double *dens) {
  int i, j, aux, counter = 0;
  BasisSet basis = sys->basis;
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

// void becke_weights(int npoint, double *points, double *weights, int natom,
//                    double *radii, double *centers, int select, int order) {
//   double nom,
//       denom;   // The nominator and the denominator in the weight definition
//   double p, s; // Used to build up the value of the switching function

//   // precompute the the alpha parameters for each atom pair
//   double alphas[(natom * (natom + 1)) / 2];
//   long offset = 0;
//   for (int iatom0 = 0; iatom0 < natom; iatom0++) {
//     for (int iatom1 = 0; iatom1 <= iatom0; iatom1++) {
//       // Heteronuclear assignment of the boundary. (Appendix in Becke's paper.)
//       double alpha = (radii[iatom0] - radii[iatom1]) /
//                      (radii[iatom0] + radii[iatom1]); // Eq. (A6)
//       alpha = alpha / (alpha * alpha - 1);            // Eq. (A5)
//       // Eq. (A3), except that we use some safe margin (0.45 instead of 0.5)
//       // to stay away from a ridiculous imbalance.
//       if (alpha > 0.45) {
//         alpha = 0.45;
//       } else if (alpha < -0.45) {
//         alpha = -0.45;
//       }
//       alphas[offset] = alpha;
//       offset += 1;
//     }
//   }

//   // precompute interatomic distances
//   double atomic_dists[(natom * (natom + 1)) / 2];
//   offset = 0;
//   for (int iatom0 = 0; iatom0 < natom; iatom0++) {
//     for (int iatom1 = 0; iatom1 <= iatom0; iatom1++) {
//       atomic_dists[offset] = dist(&centers[3 * iatom0], &centers[3 * iatom1]);
//       offset += 1;
//     }
//   }

//   // actual computations of Becke weights
//   for (int ipoint = npoint - 1; ipoint >= 0; ipoint--) {
//     nom = 0;
//     denom = 0;
//     for (int iatom0 = 0; iatom0 < natom; iatom0++) {
//       p = 1;
//       for (int iatom1 = 0; iatom1 < natom; iatom1++) {
//         if (iatom0 == iatom1)
//           continue;

//         // compute offset for alpha and interatomic distance
//         if (iatom0 < iatom1) {
//           offset = (iatom1 * (iatom1 + 1)) / 2 + iatom0;
//         } else {
//           offset = (iatom0 * (iatom0 + 1)) / 2 + iatom1;
//         }

//         // Diatomic switching function
//         s = (dist(points, &centers[3 * iatom0]) -
//              dist(points, &centers[3 * iatom1])) /
//             atomic_dists[offset]; // Eq. (11)
//         s = s +
//             alphas[offset] * (1 - 2 * (iatom0 < iatom1)) *
//                 (1 - s * s); // Eq. (A2)

//         for (int k = 1; k <= order; k++) { // Eq. (19) and (20)
//           s = 0.5 * s * (3 - s * s);
//         }
//         s = 0.5 * (1 - s); // Eq. (18)

//         p *= s; // Eq. (13)
// #ifdef DEBUG
//         printf("iatom0=%i  iatom1=%i s=%f p=%f\n", iatom0, iatom1, s, p);
// #endif
//       }

//       if (iatom0 == select)
//         nom = p;
//       denom += p; // Eq. (22)
//     }
// #ifdef DEBUG
//     printf("nom=%f  denom=%f\n", nom, denom);
// #endif

//     // Weight function at this grid point:
//     *weights *= nom / denom; // Eq. (22)

//     // go to next point
//     points += 3;
//     weights++;
//   }
// }
