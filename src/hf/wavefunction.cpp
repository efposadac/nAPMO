/*
file: wavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co
*/

#include "wavefunction.h"

void wavefunction_guess_hcore(WaveFunction *psi) {

  int ndim = psi->ndim;

  MMap S(psi->S, ndim, ndim);
  MMap H(psi->H, ndim, ndim);
  MMap C(psi->C, ndim, ndim);
  MMap D(psi->D, ndim, ndim);
  VMap O(psi->O, ndim);

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);

  O = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tHcore Matrix: " << "occupation: "<< psi->occupation<<" "<<psi->eta<<"\n";
  // std::cout << H << std::endl;

  // compute density, D = C(occ) . C(occ)T
  auto C_occ = C.leftCols(psi->occupation);
  D = C_occ * C_occ.transpose();
  D *= psi->eta;

  // std::cout.precision(20);
  // std::cout << "\n\tDensity Matrix: "<<D.sum()<<"\n";
  // std::cout << D << std::endl;
}

void wavefunction_transformation_matrix(WaveFunction *psi) {
  int ndim = psi->ndim;
  MMap X(psi->X, ndim, ndim);
  VMap O(psi->O, ndim);

  // Eigen::EigenSolver<Matrix> eig_solver(S);

  // auto eig_val = eig_solver.eigenvalues().real();
  Matrix eig_vec = X;

  using namespace std;

  // cout << O.sum() << " " << X.sum() << endl;

  for (int i = 0; i < ndim; ++i)
    for (int j = 0; j < ndim; ++j) {
      X(i, j) /= sqrt(O(j));
    }

  // cout << X.sum() << " "<< eig_vec.sum()<< endl;

  X *= eig_vec.transpose();

}


void wavefunction_iterate(WaveFunction *psi) {

  int ndim = psi->ndim;
  double xc_energy = psi->xc_energy;

  MMap S(psi->S, ndim, ndim);
  MMap C(psi->C, ndim, ndim);
  MMap H(psi->H, ndim, ndim);
  MMap D(psi->D, ndim, ndim);
  MMap L(psi->L, ndim, ndim);
  MMap G(psi->G, ndim, ndim);
  MMap J(psi->J, ndim, ndim);
  MMap F(psi->F, ndim, ndim);
  VMap O(psi->O, ndim);

  L = D;

  // std::cout << "\n\tFock Matrix:\n";
  // std::cout << F << std::endl;

  // std::cout << "\n\tS Matrix:\n";
  // std::cout << S << std::endl;

  // std::cout << "\n\tDensity Matrix:\n";
  // std::cout << D << std::endl;

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);

  O = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tEigenvectors:\n";
  // std::cout << O << std::endl;

  // std::cout << "\n\tC Matrix:\n";
  // std::cout << C << std::endl;

  // compute density, D = C(occ) . C(occ)T
  auto C_occ = C.leftCols(psi->occupation);
  D = C_occ * C_occ.transpose();
  D *= psi->eta;

  // std::cout << "\n\tDensity Matrix:\n";
  // std::cout << D << std::endl;

  // compute HF energy
  auto ehf = D.cwiseProduct(H + (0.5 * G) + J).sum() + xc_energy ;
  psi->energy = ehf;
  psi->rmsd = (D - L).norm();
}

void wavefunction_compute_coefficients(WaveFunction *psi) {

  int ndim = psi->ndim;
  MMap S(psi->S, ndim, ndim);
  MMap C(psi->C, ndim, ndim);
  MMap F(psi->F, ndim, ndim);

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
  auto eps = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tC Matrix:\n";
  // std::cout << C << std::endl;
}

void wavefunction_compute_density(WaveFunction *psi) {

  int ndim = psi->ndim;
  MMap C(psi->C, ndim, ndim);
  MMap D(psi->D, ndim, ndim);
  MMap L(psi->L, ndim, ndim);

  L = D;

  // compute density, D = C(occ) . C(occ)T
  auto C_occ = C.leftCols(psi->occupation);
  D = C_occ * C_occ.transpose();
  D *= psi->eta;

  // std::cout << "\n\tDensity Matrix:\n";
  // std::cout << D << std::endl;

  psi->rmsd = (D - L).norm();
}

void wavefunction_compute_energy(WaveFunction *psi) {

  int ndim = psi->ndim;
  double xc_energy = psi->xc_energy;

  MMap H(psi->H, ndim, ndim);
  MMap D(psi->D, ndim, ndim);
  MMap G(psi->G, ndim, ndim);
  MMap J(psi->J, ndim, ndim);

  // compute HF energy
  auto ehf = D.cwiseProduct(H + (0.5 * G) + J).sum() + xc_energy;
  psi->energy = ehf;
}

void wavefunction_compute_2body_matrix(WaveFunction *psi,
                                       std::vector<QuartetBuffer> *ints) {

  using napmo::nthreads;

  set_nthreads();

  int ndim = psi->ndim;
  double factor = psi->x_factor;

  MMap D(psi->D, ndim, ndim);
  MMap G(psi->G, ndim, ndim);

  G.setZero();

  // FELIX: add a conditional for HF or hybrid functionals
  // printf("exchange factor %f \n", factor);

  std::vector<Matrix> GB(nthreads, Matrix::Zero(ndim, ndim));

  auto lambda = [&](unsigned int thread_id) {
    auto &g = GB[thread_id];

    for (unsigned int i = 0; i < ints->at(thread_id).p.size(); i++) {

      auto p = ints->at(thread_id).s[i];
      auto q = ints->at(thread_id).r[i];
      auto r = ints->at(thread_id).q[i];
      auto s = ints->at(thread_id).p[i];
      auto val = ints->at(thread_id).val[i];

      // std::cout << p << " " << q << " " << r << " " << s << " " << val << std::endl;

      auto coulomb = D(r, s) * val;

      // Adds coulomb operator contributions

      if (p == r && q == s) {
        g(p, q) = g(p, q) + coulomb;
        if (r != s) {
          g(p, q) = g(p, q) + coulomb;
        }
      } else {
        g(p, q) = g(p, q) + coulomb;
        if (r != s) {
          g(p, q) = g(p, q) + coulomb;
        }
        coulomb = D(p, q) * val;
        g(r, s) = g(r, s) + coulomb;
        if (p != q) {
          g(r, s) = g(r, s) + coulomb;
        }
      }

      // Adds exchange operator contributions
      auto exchange = 0.0;
      if (factor != 0.0) {
        if (r != s) {
          exchange = D(q, s) * val * factor;
          g(p, r) += exchange;

          if (p == r && q != s) {
            g(p, r) += exchange;
          }
        }
        if (p != q) {
          exchange = D(p, r) * val * factor;
          if (q > s) {
            g(s, q) += exchange;
          } else {
            g(q, s) += exchange;

            if (q == s && p != r) {
              g(q, s) += exchange;
            }
          }
          if (r != s) {
            exchange = D(p, s) * val * factor;
            if (q <= r) {
              g(q, r) += exchange;
              if (q == r) {
                g(q, r) += exchange;
              }
            } else {
              g(r, q) += exchange;
              if (p == r && s == q)
                continue;
            }
          }
        }

        exchange = D(q, r) * val * factor;
        g(p, s) += exchange;
      }
    }
  }; // end lambda

  napmo::parallel_do(lambda);

  // accumulate contributions from all threads
  for (size_t i = 1; i < nthreads; ++i) {
    GB[0] += GB[i];
  }

  // Make symmetric
  GB[0] += GB[0].triangularView<Eigen::StrictlyLower>().transpose();
  G = GB[0].triangularView<Eigen::Upper>();
  G += G.triangularView<Eigen::StrictlyUpper>().transpose();
}

