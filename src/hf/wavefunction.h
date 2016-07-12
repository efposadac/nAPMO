/*
file: wavefunction.h
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

// #include <stdlib.h>
// #include <iostream>

// Eigen matrix algebra library
// #include <Eigen/Dense>
// #include <Eigen/Eigenvalues>

#include "../ints/ints.h"

struct wf {
  double *S; // Overlap matrix
  double *T; // Kinetic matrix
  double *V; // Nuclear matrix
  double *H; // One-particle Hamiltonian
  double *C; // Coefficients Matrix
  double *D; // Density matrix
  double *L; // Last Density matrix
  double *G; // two-particle matrix
  double *J; // Coupling matrix
  double *F; // Fock matrix
  // Convenience Variables
  int nbasis; // order of matrices
  int occupation;
  double eta;    // Constant of coupling
  double kappa;  // Constant of coupling
  double energy; // HF Energy
  double rmsd;   // Density root-mean-square deviation
};
typedef struct wf WaveFunction;

/*
Type definitions
*/
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix; // import dense, dynamically sized Matrix type from Eigen;
            // this is a matrix with row-major storage
            // (http://en.wikipedia.org/wiki/Row-major_order)
typedef Eigen::Map<Matrix> Map;

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

#ifdef __cplusplus
extern "C" {
#endif

void wavefunction_guess_hcore(WaveFunction *psi);

void wavefunction_compute_density(WaveFunction *psi);

void wavefunction_iterate(WaveFunction *psi);

void wavefunction_compute_energy(WaveFunction *psi);

void wavefunction_compute_2body_matrix(WaveFunction *psi,
                                       std::vector<QuartetBuffer> *ints);

#ifdef __cplusplus
}
#endif

#endif
