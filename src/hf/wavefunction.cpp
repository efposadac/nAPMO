/*
file: wavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"

namespace napmo {

unsigned int nthreads;

// / fires off \c nthreads instances of lambda in parallel
template <typename Lambda> void parallel_do(Lambda &lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else // use C++11 threads
  std::vector<std::thread> threads;
  for (unsigned int thread_id = 0; thread_id != napmo::nthreads; ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  } // threads_id
  for (unsigned int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}
} // end namespace

__inline void set_nthreads() {
  {
    using napmo::nthreads;
    auto nthreads_cstr = getenv("OMP_NUM_THREADS");
    nthreads = 1;
    if (nthreads_cstr && strcmp(nthreads_cstr, "")) {
      std::istringstream iss(nthreads_cstr);
      iss >> nthreads;
      if (nthreads > 1 << 16 || nthreads <= 0)
        nthreads = 1;
    }

#if defined(_OPENMP)
    omp_set_num_threads(nthreads);
#endif
    //     std::cout << "Will scale over " << nthreads
    // #if defined(_OPENMP)
    //               << " OpenMP"
    // #else
    //               << " C++11"
    // #endif
    //               << " threads" << std::endl;
  }
}

void wavefunction_guess_hcore(WaveFunction *psi) {

  int nbasis = psi->nbasis;

  Map S(psi->S, nbasis, nbasis);
  Map H(psi->H, nbasis, nbasis);
  Map C(psi->C, nbasis, nbasis);
  Map D(psi->D, nbasis, nbasis);

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
  auto eps = gen_eig_solver.eigenvalues();

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
  Map S(psi->S, nbasis, nbasis);
  Map C(psi->C, nbasis, nbasis);
  Map H(psi->H, nbasis, nbasis);
  Map D(psi->D, nbasis, nbasis);
  Map L(psi->L, nbasis, nbasis);
  Map G(psi->G, nbasis, nbasis);
  Map F(psi->F, nbasis, nbasis);

  L = D;

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
  auto eps = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

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
  Map S(psi->S, nbasis, nbasis);
  Map C(psi->C, nbasis, nbasis);
  Map F(psi->F, nbasis, nbasis);

  Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);
  auto eps = gen_eig_solver.eigenvalues();
  C = gen_eig_solver.eigenvectors();

  // std::cout << "\n\tC Matrix:\n";
  // std::cout << C << std::endl;
}

void wavefunction_compute_density(WaveFunction *psi) {

  int nbasis = psi->nbasis;
  Map C(psi->C, nbasis, nbasis);
  Map D(psi->D, nbasis, nbasis);
  Map L(psi->L, nbasis, nbasis);

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
  Map H(psi->H, nbasis, nbasis);
  Map D(psi->D, nbasis, nbasis);
  Map G(psi->G, nbasis, nbasis);
  Map J(psi->J, nbasis, nbasis);
  
  // compute HF energy
  auto ehf = D.cwiseProduct(H + (0.5 * G) + J).sum();
  psi->energy = ehf;
}

void wavefunction_compute_2body_matrix(WaveFunction *psi,
                                       std::vector<QuartetBuffer> *ints) {

  using napmo::nthreads;

  set_nthreads();

  int nbasis = psi->nbasis;
  Map D(psi->D, nbasis, nbasis);
  Map G(psi->G, nbasis, nbasis);

  G.setZero();

  auto factor = psi->kappa / psi->eta;

  std::vector<Matrix> GB(nthreads, Matrix::Zero(nbasis, nbasis));

  auto lambda = [&](unsigned int thread_id) {

    auto &g = GB[thread_id];

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
        g(p, r) = g(p, r) + exchange;
        if (p == r && q != s) {
          g(p, r) = g(p, r) + exchange;
        }
      }
      if (p != q) {
        exchange = D(p, r) * val * factor;
        if (q > s) {
          g(s, q) = g(s, q) + exchange;
        } else {
          g(q, s) = g(q, s) + exchange;
          if (q == s && p != r) {
            g(q, s) = g(q, s) + exchange;
          }
        }
        if (r != s) {
          exchange = D(p, s) * val * factor;
          if (q <= r) {
            g(q, r) = g(q, r) + exchange;
            if (q == r) {
              g(q, r) = g(q, r) + exchange;
            }
          } else {
            g(r, q) = g(r, q) + exchange;
            if (p == r && s == q)
              continue;
          }
        }
      }

      exchange = D(q, r) * val * factor;
      g(p, s) = g(p, s) + exchange;
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