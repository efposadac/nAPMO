/*
file: nwavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"

#define INDEX(i, j)                                                            \
  ((i > j) ? (((i) * ((i) + 1) / 2) + (j)) : (((j) * ((j) + 1) / 2) + (i)))

void nwavefunction_compute_2body_matrix(WaveFunction *psi, BeckeGrid *grid,
                                        double *phi, double *J, double *K) {

  unsigned int ndim = psi->ndim;
  unsigned int size = ndim * (ndim + 1) / 2;

  // Precompute Psi
  MMap F(phi, ndim, grid->get_size());
  Matrix Psi(size, grid->get_size());

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
    }
  }

  // Compute coulomb
  A1DMap C(J, grid->get_size());
  Matrix coulomb(ndim, ndim);
  coulomb.setZero();

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= C;

      coulomb(i, j) += grid->integrate(buff);
    }
  }

  coulomb += coulomb.triangularView<Eigen::StrictlyUpper>().transpose();

  // Compute exchange
  MMap D(psi->D, ndim, ndim);
  MMap E(K, size, grid->get_size());
  Matrix exchange(ndim, ndim);
  exchange.setZero();

  double factor = psi->kappa / psi->eta;

  for (unsigned int u = 0; u < ndim; ++u) {
    for (unsigned int v = u; v < ndim; ++v) {
      for (unsigned int l = 0; l < ndim; ++l) {
        for (unsigned int s = 0; s < ndim; ++s) {

          Array1D buff =
              Psi.row(INDEX(u, l)).array() * E.row(INDEX(s, v)).array();

          exchange(u, v) += grid->integrate(buff) * D(l, s);
        }
      }
    }
  }

  // MMap E(K, ndim, grid->get_size());
  // Matrix exchange(ndim, ndim);
  // exchange.setZero();

  // double factor = psi->kappa / psi->eta;
  // for (unsigned int i = 0; i < ndim; ++i) {
  //   for (unsigned int j = i; j < ndim; ++j) {

  //     Array1D buff = F.row(i).array() * E.row(j).array();

  //     exchange(i, j) += grid->integrate(buff);
  //   }
  // }

  exchange += exchange.triangularView<Eigen::StrictlyUpper>().transpose();
  // exchange *= factor * psi->eta;
  exchange *= factor;
  
  // std::cout<<exchange<<std::endl;

  MMap G(psi->G, ndim, ndim);

  G = coulomb + exchange;
}

void nwavefunction_compute_density_from_dm(BasisSet *basis, BeckeGrid *grid,
                                           double *dm, double *output,
                                           double epsilon, double *dmmaxrow) {
  int nbasis = basis->get_nbasis();

  for (int point = 0; point < grid->get_size(); ++point) {

    auto aux = basis->compute(&grid->get_points()[point * 3]);

    if (epsilon > 0) {

      // compute the maximum basis function
      double absmax_basis = 0.0;
      for (int i = 0; i < nbasis; i++) {
        absmax_basis = std::max(fabs(aux[i]), absmax_basis);
      }

      // upper estimate of the density
      double rho_upper = 0.0;
      for (int i = 0; i < nbasis; i++) {
        rho_upper += fabs(aux[i]) * dmmaxrow[i];
      }
      rho_upper *= nbasis * absmax_basis;

      // if the upper bound is too low, do not compute density.
      if (rho_upper < epsilon)
        return;

      // modify epsilon to avoid recomputation
      epsilon /= absmax_basis * nbasis * nbasis;
    }

    // Loop over all basis functions and add significant contributions
    double rho = 0.0;
    for (int i = 0; i < nbasis; i++) {

      // if the contribution of this loop is smaller than
      // epsilon/nbasis, skipt
      // it.
      if (epsilon > 0) {
        if (fabs(aux[i]) * dmmaxrow[i] < epsilon)
          continue;
      }

      double tmp = 0;
      for (int j = i - 1; j >= 0; j--) {
        tmp += aux[j] * dm[i * nbasis + j];
      }
      rho += (2.0 * tmp + dm[i * (nbasis + 1)] * aux[i]) * aux[i];
    }

    output[point] = rho;
  }
}

////////////////////////////TESTING///////////////////////////

double nwavefunction_energy(double a, double b, double c1, Matrix E) {
  return (pow(a, 2) * E(0, 0) + pow(b, 2) * E(1, 1) -
          sqrt(pow(a, 4) * pow(E(0, 0), 2) + pow(b, 4) * pow(E(1, 1), 2) +
               2 * pow(a, 2) * pow(b, 2) *
                   (2 * pow(E(0, 1), 2) - E(0, 0) * E(1, 1)))) /
         2.;
}

Vector nwavefunction_energy_grads(double a, double b, double c1, Matrix E,
                                  Matrix s) {
  Vector result(3);

  result[0] =
      (2 * a * E(0, 0) -
       (2 * a * (pow(a, 2) * pow(E(0, 0), 2) +
                 pow(b, 2) * (2 * pow(E(0, 1), 2) - E(0, 0) * E(1, 1)))) /
           sqrt(pow(a, 4) * pow(E(0, 0), 2) + pow(b, 4) * pow(E(1, 1), 2) +
                2 * pow(a, 2) * pow(b, 2) *
                    (2 * pow(E(0, 1), 2) - E(0, 0) * E(1, 1)))) /
          2. -
      2 * c1 * (a * s(0, 0) + b * s(0, 1));

  result[1] =
      (2 * b * E(1, 1) -
       (2 * b * (pow(b, 2) * pow(E(1, 1), 2) +
                 pow(a, 2) * (2 * pow(E(0, 1), 2) - E(0, 0) * E(1, 1)))) /
           sqrt(pow(a, 4) * pow(E(0, 0), 2) + pow(b, 4) * pow(E(1, 1), 2) +
                2 * pow(a, 2) * pow(b, 2) *
                    (2 * pow(E(0, 1), 2) - E(0, 0) * E(1, 1)))) /
          2. -
      2 * c1 * (a * s(0, 1) + b * s(1, 1));

  result[2] = 1 - a * (a * s(0, 0) + 2 * b * s(0, 1)) - pow(b, 2) * s(1, 1);

  return result;
}

Vector nwavefunction_orbital_optimizer(Vector a, int n_steps, double gamma,
                                       Matrix E, Matrix s) {
  Vector ap = a, b;
  int counter = 0;
  double old_v, new_v;

  while (counter < n_steps) {
    old_v = nwavefunction_energy(ap[0], ap[1], ap[2], E);
    b = ap - gamma * nwavefunction_energy_grads(ap[0], ap[1], ap[2], E, s);
    ap = b;
    new_v = nwavefunction_energy(ap[0], ap[1], ap[2], E);

    printf("%d %e %e\n", counter, nwavefunction_energy(ap[0], ap[1], ap[2], E),
           abs(old_v - new_v));
    counter++;
    // if(abs(old_v-new_v)<1.0e-10)
    //   break;
  }
  return ap;
}

void nwavefunction_optimize(double *H, double *S, double *G, double *res,
                            int ndim) {
  MMap E(H, ndim, ndim);
  MMap s(S, ndim, ndim);
  VMap R(res, 3);

  Vector guess(3);
  guess << G[0], G[1], 0.95;

  std::cout << E << std::endl;
  std::cout << std::endl;

  std::cout << s << std::endl;
  std::cout << std::endl;

  std::cout << guess << std::endl;

  R = nwavefunction_orbital_optimizer(guess, 500, 0.00001, E, s);
  std::cout << R << std::endl;
}