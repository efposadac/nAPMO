/*
file: nwavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"

void nwavefunction_compute_2body_matrix(WaveFunction *psi, double *func,
                                        int size) {
  auto nbasis = psi->nbasis;

  Matrix coulomb(nbasis, nbasis);
  MMap F(func, nbasis, size);

  for (int i = 0; i < nbasis; ++i) {
    for (int j = 0; j < nbasis; ++j) {
      // coulomb[i, j] += self._grid.integrate(
      //         self.psi[i, :] * self.psi[j, :] * self.JO)
    }
  }

  // factor = self._kappa / self._eta
  // exchange = np.zeros([self.nbasis, self.nbasis])

  // for s in range(self.nbasis):
  //   for v in range(self.nbasis):
  //       for u in range(self.nbasis):
  //           for l in range(self.nbasis):
  //               exchange[
  //                   u, v] += self._grid.integrate(
  //                   self.psi[u, :] * self.psi[l, :] *
  //                   self.KO[index2(s, v)]) * self.D[l, s]

  // exchange *= factor

  // self.G[:] = coulomb + exchange
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

      // if the contribution of this loop is smaller than epsilon/nbasis, skipt
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