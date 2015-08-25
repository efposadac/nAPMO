/*file: test_gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "test.h"

void test_gauss_chebyshev() {
  int status = 1;
  int i, n;
  double eps = 1e-8;
  double *abscissas;
  double *weights;

  double _r[] = {1.4179599213, 0.3803914706, 0.0645385211, 0.0039138993, 0.0000201360};
  double _w[] = {1.7557935946, 0.5596866610, 0.1396263402, 0.0142258774, 0.0001573928};

  n = 5;

  abscissas = (double *)malloc(n * sizeof(double));
  weights = (double *)malloc(n * sizeof(double));

  gaussChebyshev(n, abscissas, weights);

  for (i = 0; i < n; ++i) {
    if (!(fabs(abscissas[i] - _r[i]) < eps)) {
      printf("Fail! expected: %12.10f got: %12.10f\n", _r[i], abscissas[i]);
      status = 0;
    }
    if (!(fabs(weights[i] - _w[i]) < eps)) {
      printf("Fail! expected: %12.10f got: %12.10f\n", _w[i], weights[i]);
      status = 0;
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

  reps = 100;

  for (n = 2; n < 1000000; n += 10000) {
    printf("%d ", n);

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
