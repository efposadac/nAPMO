/*
file: wavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"

void wavefunction_guess_hcore(WaveFunction *psi) {

  int nbasis = psi->nbasis;

  MMap S(psi->S, nbasis, nbasis);
  MMap H(psi->H, nbasis, nbasis);
  MMap C(psi->C, nbasis, nbasis);
  MMap D(psi->D, nbasis, nbasis);
  VMap O(psi->O, nbasis);

  // std::cout<<nbasis<<psi->occupation<<"\n";

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);

  O = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tC Matrix: " << "occupation: "<< psi->occupation<<"\n";
  // std::cout << C << std::endl;

  // compute density, D = C(occ) . C(occ)T
  auto C_occ = C.leftCols(psi->occupation);
  D = C_occ * C_occ.transpose();
  D *= psi->eta;

  // std::cout << "\n\tDensity Matrix:\n";
  // std::cout << D << std::endl;
}

void wavefunction_iterate(WaveFunction *psi) {

  int nbasis = psi->nbasis;
  MMap S(psi->S, nbasis, nbasis);
  MMap C(psi->C, nbasis, nbasis);
  MMap H(psi->H, nbasis, nbasis);
  MMap D(psi->D, nbasis, nbasis);
  MMap L(psi->L, nbasis, nbasis);
  MMap G(psi->G, nbasis, nbasis);
  MMap F(psi->F, nbasis, nbasis);
  VMap O(psi->O, nbasis);

  L = D;

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
  auto ehf = D.cwiseProduct(H + (0.5 * G)).sum();
  psi->energy = ehf;
  psi->rmsd = (D - L).norm();
}

void wavefunction_compute_coefficients(WaveFunction *psi) {

  int nbasis = psi->nbasis;
  MMap S(psi->S, nbasis, nbasis);
  MMap C(psi->C, nbasis, nbasis);
  MMap F(psi->F, nbasis, nbasis);

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
  auto eps = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tC Matrix:\n";
  // std::cout << C << std::endl;
}

void wavefunction_compute_density(WaveFunction *psi) {

  int nbasis = psi->nbasis;
  MMap C(psi->C, nbasis, nbasis);
  MMap D(psi->D, nbasis, nbasis);
  MMap L(psi->L, nbasis, nbasis);

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

  int nbasis = psi->nbasis;
  MMap H(psi->H, nbasis, nbasis);
  MMap D(psi->D, nbasis, nbasis);
  MMap G(psi->G, nbasis, nbasis);
  MMap J(psi->J, nbasis, nbasis);

  // compute HF energy
  auto ehf = D.cwiseProduct(H + (0.5 * G) + J).sum();
  psi->energy = ehf;
}

void wavefunction_compute_2body_matrix(WaveFunction *psi,
                                       std::vector<QuartetBuffer> *ints) {

  using napmo::nthreads;

  set_nthreads();

  int nbasis = psi->nbasis;
  MMap D(psi->D, nbasis, nbasis);
  MMap G(psi->G, nbasis, nbasis);

  G.setZero();

  auto factor = psi->kappa / psi->eta;

  std::vector<Matrix> GB(nthreads, Matrix::Zero(nbasis, nbasis));
  std::vector<Matrix> TEST(nthreads, Matrix::Zero(nbasis, nbasis));

  auto lambda = [&](unsigned int thread_id) {

    auto &g = GB[thread_id];
    auto &t = TEST[thread_id];

    t.setZero();

    for (unsigned int i = 0; i < ints->at(thread_id).p.size(); i++) {

      auto p = ints->at(thread_id).s[i];
      auto q = ints->at(thread_id).r[i];
      auto r = ints->at(thread_id).q[i];
      auto s = ints->at(thread_id).p[i];
      auto val = ints->at(thread_id).val[i];

      // cout << p << " " << q << " " << r << " " << s << " " << val << endl;

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