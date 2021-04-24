/*file: auxiliary_basis.cpp
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#include "auxiliary_basis.h"

AuxiliaryBasis::AuxiliaryBasis(BeckeGrid *grid) {
  molgrid = grid;
  ncenter = molgrid->get_ncenter();
  abldep = molgrid->get_abldep(); // Threshold for linear dependency
  lcenter.reserve(ncenter);

  ablmax = molgrid->get_ablmax();

  printf("L max for Auxiliary Basis: %d\n", ablmax);
  printf("Linear Dependency Tolerance: %12.9E\n", abldep);
  printf("Number of Centers: %d\n", ncenter);

  max_rad = 0; // largest radial quadrature
  max_atm = 0; // largest atomic grid

  for (int center = 0; center < ncenter; ++center) {
    max_rad = std::max(molgrid->atgrid[center].rad_grid->get_size(), max_rad);
    max_atm = std::max(molgrid->atgrid[center].get_size(), max_atm);
  }

  compute_aobasis();
}

void AuxiliaryBasis::compute_aobasis() {
  // Determine exponents of the Slater-type auxiliary basis functions
  // s-, p-, d-, ..., AG-type real Cartesian Slater-spherical harmonics

  std::vector<int> nbas(ncenter);
  std::vector<Array2D> aoexp(ncenter, Vector(max_rad));
  std::vector<Matrix> C(ncenter, Matrix(max_rad, max_rad));

  // Remove linear dependencies for s-type for each center
  int nbastotal = 0;
  double log_2 = std::log(2.0);
  double sqrt_2 = std::sqrt(2.0);

  for (int ip = 0; ip < ncenter; ++ip) {
    AtomicGrid *atgrid = &molgrid->atgrid[ip];

    unsigned int asize = atgrid->get_size();
    A2DMap apoints(atgrid->get_points(), 3, asize);
    A1DMap aweights(atgrid->get_weights(), asize);

    unsigned int rsize = atgrid->rad_grid->get_size();
    A1DMap rpoints(atgrid->rad_grid->get_points(), rsize);

    // Determine the exponents of slater auxiliary functions
    // Gauss-Chevyshev radial-grid based
    aoexp[ip] = log_2 / rpoints;

    // Build overlap matrix
    Vector RR(asize);
    // auto aux = apoints - atgrid->get_origin();
    RR = Eigen::sqrt((apoints * apoints).colwise().sum());
    Matrix A(rsize, rsize);
    A.setZero();
    for (unsigned int i = 0; i < asize; ++i) {
      Vector tmp = Eigen::exp(-aoexp[ip] * RR(i));
      A += tmp * (tmp.transpose() * aweights(i));
    }

    // Diagonalize overlap && sort eigenvalues
    Matrix eigenvectors(rsize, rsize);
    Vector eigenvalues(rsize);
    A.eigenVectorsVec(eigenvectors, eigenvalues);
    gram_schmidt(&eigenvectors);

    // std::cout<<std::fixed<<std::setprecision(12)<<eigenvectors<<std::endl;

    // canonical orthogonalization
    nbas[ip] = 0;
    for (unsigned int i = 0; i < rsize; ++i) {
      if (eigenvalues(i) > abldep) {
        C[ip].col(nbas[ip]) = eigenvectors.col(i) / std::sqrt(eigenvalues(i));
        nbas[ip]++;
      }
    }

    printf("Number of exponents for particle %d: %d\n", ip, nbas[ip]);
    C[ip].conservativeResize(rsize, nbas[ip]);

    nbastotal += nbas[ip];
  }

  // Construct real-space AO
  int norb = 0;
  double r;
  std::complex<double> ylm1;
  std::complex<double> ylm2;
  std::vector<Matrix> aotemp(
      ncenter, Matrix(max_atm, int(nbastotal * std::pow(ablmax + 1, 2))));

  for (int ip = 0; ip < ncenter; ++ip) {
    AtomicGrid *atgrid_i = &molgrid->atgrid[ip];
    A1DMap origin_i(atgrid_i->get_origin(), 3);

    // S-type AO
    for (int jp = 0; jp < ncenter; ++jp) {
      AtomicGrid *atgrid_j = &molgrid->atgrid[jp];

      unsigned int size_j = atgrid_j->get_size();
      A2DMap points_j(atgrid_j->get_points(), 3, size_j);

      Vector RR_j(size_j);
      // auto aux = points_j - origin_i;
      RR_j = Eigen::sqrt((points_j * points_j).colwise().sum());

      for (unsigned int i = 0; i < size_j; ++i) {
        Vector tmp = Eigen::exp(-aoexp[ip] * RR_j(i));
        aotemp[jp].block(i, norb, 1, nbas[ip]) = tmp.transpose() * C[ip];
      }
      norb += nbas[ip];
    }

    // Populate p-type and higher AO
    if (ablmax > 0) {
      for (int l = 1; l <= ablmax; ++l) {
        for (int jp = 0; jp < ncenter; ++jp) {
          AtomicGrid *atgrid_j = &molgrid->atgrid[jp];

          unsigned int size_j = atgrid_j->get_size();
          A2DMap points_j(atgrid_j->get_points(), 3, size_j);

          for (unsigned int k = 0; k < size_j; ++k) {
            rylm(points_j.col(k), origin_i, l, 0, &ylm1, &r);

            double factor = std::real(ylm1);
            aotemp[jp].block(k, norb, 1, nbas[ip]) =
                aotemp[jp].block(k, 0, 1, nbas[ip]) * factor;
          }
          norb += nbas[ip];
        }
        for (int m = 1; m <= l; ++m) {
          double sign = (m % 2) ? -1.0 : 1.0;
          for (int jp = 0; jp < ncenter; ++jp) {
            AtomicGrid *atgrid_j = &molgrid->atgrid[jp];

            unsigned int size_j = atgrid_j->get_size();
            A2DMap points_j(atgrid_j->get_points(), 3, size_j);

            for (unsigned int k = 0; k < size_j; ++k) {
              rylm(points_j.col(k), origin_i, l, +m, &ylm1, &r);
              rylm(points_j.col(k), origin_i, l, -m, &ylm2, &r);

              double factor1 = std::real(ylm1 + sign * ylm2) / sqrt_2;
              double factor2 = std::imag(ylm1 - sign * ylm2) / sqrt_2;

              aotemp[jp].block(k, norb, 1, nbas[ip]) =
                  aotemp[jp].block(k, 0, 1, nbas[ip]) * factor1;

              aotemp[jp].block(k, norb + nbas[ip], 1, nbas[ip]) =
                  aotemp[jp].block(k, 0, 1, nbas[ip]) * factor2;
            }
            norb += nbas[ip] * 2;
          }
        }
      }
    }
  }

  // Build overlap matrix
  Matrix AO(norb, norb);
  AO.setZero();
  for (int ip = 0; ip < ncenter; ++ip) {
    AtomicGrid *atgrid = &molgrid->atgrid[ip];
    auto asize = atgrid->get_size();
    auto aweights = atgrid->get_weights();

    for (unsigned int k = 0; k < asize; ++k) {
      AO += aotemp[ip].row(k).transpose() * aotemp[ip].row(k) * aweights[k];
    }
  }

  // Diagonalize overlap && sort the eigenvalues
  Vector eigenvalues(norb);
  Matrix eigenvectors(norb, norb);
  AO.eigenVectorsVec(eigenvectors, eigenvalues);
  gram_schmidt(&eigenvectors);

  // canonical orthogonalization
  Matrix aux(norb, norb);
  int n = 0;
  for (int i = 0; i < norb; ++i) {
    // remove linear dependencies
    if (eigenvalues(i) > abldep) {
      for (int j = 0; j < norb; ++j) {
        aux(j, n) = eigenvectors(j, i) / sqrt(eigenvalues(i));
      }
      n++;
    }
  }

  printf("Number of Linear Dependencies (LMAX: %d): %d\n", ablmax, norb - n);
  printf("Number of Auxiliary Functions (LMAX: %d): %d\n", ablmax, n);

  nao = n;
  aobasis = new double[ncenter * max_atm * nao];

  for (int ip = 0; ip < ncenter; ++ip) {
    AtomicGrid *atgrid = &molgrid->atgrid[ip];
    auto asize = atgrid->get_size();
    int ii = ip * asize * n;
    for (unsigned int j = 0; j < asize; ++j) {
      for (int i = 0; i < n; ++i) {
        int idx = ii + (j * n + i);
        aobasis[idx] = 0.0;
        for (int k = 0; k < norb; ++k) {
          aobasis[idx] += aotemp[ip](j, k) * aux(k, i);
        }
      }
    }
  }

  // Print Aux basis
  // for (int ia = 0; ia < ncenter; ++ia) {
  //   AtomicGrid *atgrid = &molgrid->atgrid[ia];
  //   auto asize = atgrid->get_size();
  //   int ii = ia * asize * n;
  //   for (unsigned int j = 0; j < asize; ++j) {
  //     for (int i = 0; i < n; ++i) {
  //       int idx = ii + (j * n + i);
  //       printf("%25.17E ", aobasis[idx]);
  //     }
  //     printf("\n");
  //   }
  // }
}

void AuxiliaryBasis::gram_schmidt(Matrix *M) {
  // Gram-Schmidt orthogonalization

  int n = M->cols();
  for (int i = 0; i < n; ++i) {
    Vector colVec = M->col(i);
    for (int j = 0; j < i; ++j) {
      Vector prevCol = M->col(j);
      double d = prevCol.dot(prevCol);
      double c = colVec.dot(prevCol);
      colVec -= c * prevCol / d;
    }
    M->col(i) = colVec.normalized();
  }
}

void AuxiliaryBasis::rylm(Array1D dest, Array1D orig, int l, int m,
                          std::complex<double> *ylm, double *r) {
  // Returns the value of radius and real spherical harmonics of the origin
  // atom at a grid point of the destination atom
  Array1D r_orig_dest(3);

  r_orig_dest = dest - orig;

  // Get Radius
  *r = std::sqrt((r_orig_dest * r_orig_dest).sum());

  // get spherical coordinates
  double t = std::acos(r_orig_dest(2) / *r);             // theta
  double p = std::atan2(r_orig_dest(1), r_orig_dest(0)); // phi

  *ylm = spherical_harmonics_complex(l, m, t, p);
}

/*
Python wrapper
*/
AuxiliaryBasis *AuxiliaryBasis_new(BeckeGrid *grid) {
  return new AuxiliaryBasis(grid);
}

void AuxiliaryBasis_del(AuxiliaryBasis *auxiliary_basis) {
  auxiliary_basis->~AuxiliaryBasis();
}

void AuxiliaryBasis_compute_aobasis(AuxiliaryBasis *auxiliary_basis) {
  auxiliary_basis->compute_aobasis();
}

double *AuxiliaryBasis_get_aobasis(AuxiliaryBasis *auxiliary_basis) {
  return auxiliary_basis->get_aobasis();
}

int AuxiliaryBasis_get_nao(AuxiliaryBasis *auxiliary_basis) {
  return auxiliary_basis->get_nao();
}

int AuxiliaryBasis_get_ncenter(AuxiliaryBasis *auxiliary_basis) {
  return auxiliary_basis->get_ncenter();
}
