/*file: omp_helper.c
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

#ifndef OMP_HELPER
#define OMP_HELPER

#include <atomic>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace napmo {

static unsigned int nthreads;

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

#ifdef _OPENMP
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

#endif