/*file: test.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "test.h"
#include "gauss_chebyshev.h"

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

  gaussChebyshev(n, abscissas, weights);

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

  // Test for performance
  n = 100000000;
  abscissas = (double *)malloc(n * sizeof(double));
  weights = (double *)malloc(n * sizeof(double));
  for (i = 0; i < 10; ++i) {
    gaussChebyshev(n, abscissas, weights);
  }
  free(abscissas);
  free(weights);
}

void test_gauss_chebyshev_perf() {
  int i, n;
  double *abscissas;
  double *weights;

  n = 100000000;
  abscissas = (double *)malloc(n * sizeof(double));
  weights = (double *)malloc(n * sizeof(double));
  for (i = 0; i < 10; ++i) {
    gaussChebyshev(n, abscissas, weights);
  }
  free(abscissas);
  free(weights);
}

int main(int argc, char const *argv[]) {
  // Results
  test_gauss_chebyshev();

  // Perf
  test_gauss_chebyshev_perf();
  return 0;
}