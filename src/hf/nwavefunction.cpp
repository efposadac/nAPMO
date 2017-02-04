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

void nwavefunction_compute_2body_matrix_atm(WaveFunction *psi, BeckeGrid *grid,
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

  exchange += exchange.triangularView<Eigen::StrictlyUpper>().transpose();
  exchange *= factor;

  // std::cout<<exchange<<std::endl;

  MMap G(psi->G, ndim, ndim);

  G = coulomb + exchange;
}

void nwavefunction_compute_2body_matrix_mol(WaveFunction *psi, BeckeGrid *grid,
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

  MMap E(K, ndim, grid->get_size());
  Matrix exchange(ndim, ndim);
  exchange.setZero();

  double factor = psi->kappa / psi->eta;
  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = F.row(i).array() * E.row(j).array();

      exchange(i, j) += grid->integrate(buff);
    }
  }

  exchange += exchange.triangularView<Eigen::StrictlyUpper>().transpose();
  exchange *= factor * psi->eta;

  // std::cout<<exchange<<std::endl;

  MMap G(psi->G, ndim, ndim);

  G = coulomb + exchange;
}

void nwavefunction_compute_coupling(WaveFunction *psi, BeckeGrid *grid,
                                    double *phi, double *other_J, double *res) {

  unsigned int ndim = psi->ndim;
  unsigned int size = ndim * (ndim + 1) / 2;

  // Precompute <ab| in the integral <ab|r12|AB> "Coupling integral"
  MMap PHI(phi, ndim, grid->get_size());
  Matrix Psi(size, grid->get_size());

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = PHI.row(i).array() * PHI.row(j).array();
    }
  }

  // Compute coupling
  A1DMap C(other_J, grid->get_size());
  MMap coupling(res, ndim, ndim);
  coupling.setZero();

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= C;

      coupling(i, j) += grid->integrate(buff);
    }
  }

  coupling += coupling.triangularView<Eigen::StrictlyUpper>().transpose();
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
