/*file: test_gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "test_gauss_chebyshev.h"

void test_gauss_chebyshev() {
  bool status = true;
  int i, n;
  double eps = 1e-8;
  double *abscissas;
  double *weights;

  double _r[] = {8.66025404e-01, 5.00000000e-01, 6.12323400e-17,
                 -5.00000000e-01, -8.66025404e-01};
  double _w[] = {0.13089969, 0.39269908, 0.52359878, 0.39269908, 0.13089969};

  n = 5;

  abscissas = (double *)malloc(n * sizeof(double));
  weights = (double *)malloc(n * sizeof(double));

#ifdef _CUDA
  gaussChebyshev_cuda(n, abscissas, weights);
#else
  gaussChebyshev(n, abscissas, weights);
#endif

  for (i = 0; i < n; ++i) {
    if (!(fabs(abscissas[i] - _r[i]) < eps)) {
      printf("Fail! expected: %12.10f got: %12.10f\n", _r[i], abscissas[i]);
      status = false;
    }
    if (!(fabs(weights[i] - _w[i]) < eps)) {
      printf("Fail! expected: %12.10f got: %12.10f\n", _w[i], weights[i]);
      status = false;
    }
  }
  printf("gauss_chebyshev... %s\n", status ? "Passed" : "Failed");

  free(abscissas);
  free(weights);
}

void test_gauss_chebyshev_perf() {
  int i, reps;
  int n;
  double *abscissas;
  double *weights;
  double secs;

#ifdef _CUDA
  double_t *g_a;
  double_t *g_w;
#endif

  reps = 100;

  for (n = 2; n < 1000000; n += 10000) {
    printf("%d ", n);

#ifdef _CUDA

    g_a = (double_t *)malloc(n * sizeof(double_t));
    g_w = (double_t *)malloc(n * sizeof(double_t));

    secs = wallclock(NULL);

    for (i = 0; i < reps; ++i) {
      gaussChebyshev_cuda(n, g_a, g_w);
    }

    free(g_a);
    free(g_w);

    printf("%g ", wallclock(&secs) / reps);
#endif

    abscissas = (double *)malloc(n * sizeof(double));
    weights = (double *)malloc(n * sizeof(double));

    secs = wallclock(NULL);

#ifdef _OMP
    omp_set_num_threads(8);
#endif

    for (i = 0; i < reps; ++i) {
      gaussChebyshev(n, abscissas, weights);
    }

#ifdef _OMP
    printf("%g ", wallclock(&secs) / reps);

    secs = wallclock(NULL);

    omp_set_num_threads(1);

    for (i = 0; i < reps; ++i) {
      gaussChebyshev(n, abscissas, weights);
    }
#endif

    printf("%g \n", wallclock(&secs) / reps);

    free(abscissas);
    free(weights);
  }
}

int main(int argc, char const *argv[]) {
  // Results
  test_gauss_chebyshev();

  // Perf
  // test_gauss_chebyshev_perf();
  return 0;
}