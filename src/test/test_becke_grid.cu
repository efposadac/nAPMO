/*file: test_becke_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "test.h"
}

void test_becke_grid_init_cuda()
{
  int status;
  int i;
  double eps = 1e-8;
  double *abscissas;
  double *weights;
  double *theta, *phi;

  /* Grid initialization*/
  Grid test_grid;
  Grid grid;

  grid.n_radial = 5;
  grid.n_angular = 110;

  grid_init(&grid);
  grid_init_cuda(&grid, &test_grid);

  /*Radial initialization check*/
  double _r[] = {8.66025404e-01, 5.00000000e-01, 6.12323400e-17, -5.00000000e-01, -8.66025404e-01};
  double _w[] = {0.13089969, 0.39269908, 0.52359878, 0.39269908, 0.13089969};

  abscissas = (double *)malloc(test_grid.n_radial * sizeof(double));
  weights = (double *)malloc(test_grid.n_radial * sizeof(double));

  cudaMemcpy(abscissas, test_grid.radial_abscissas, test_grid.n_radial * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUERR
  cudaMemcpy(weights, test_grid.radial_weights, test_grid.n_radial * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUERR

  status = 1;
  for (i = 0; i < test_grid.n_radial; ++i)
  {
    if (!(fabs(abscissas[i] - _r[i]) < eps))
    {
      printf("Fail! expected: %12.10f got: %12.10f\n", _r[i], abscissas[i]);
      status = 0;
    }
    if (!(fabs(weights[i] - _w[i]) < eps))
    {
      printf("Fail! expected: %12.10f got: %12.10f\n", _w[i], weights[i]);
      status = 0;
    }
  }
  printf("gauss_chebyshev... %s\n", status ? "Passed" : "Failed");

  free(abscissas);
  free(weights);

  /*Angular initialization check*/
  double *_at = (double *)malloc(test_grid.n_angular * sizeof(double));
  double *_ap = (double *)malloc(test_grid.n_angular * sizeof(double));
  double *_aw = (double *)malloc(test_grid.n_angular * sizeof(double));

  lebedev(test_grid.n_angular, _at, _ap, _aw);

  theta = (double *)malloc(test_grid.n_angular * sizeof(double));
  phi = (double *)malloc(test_grid.n_angular * sizeof(double));
  weights = (double *)malloc(test_grid.n_angular * sizeof(double));

  cudaMemcpy(theta, test_grid.angular_theta, test_grid.n_angular * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUERR

  cudaMemcpy(phi, test_grid.angular_phi, test_grid.n_angular * sizeof(double), cudaMemcpyDeviceToHost);
  CUERR

  cudaMemcpy(weights, test_grid.angular_weights, test_grid.n_angular * sizeof(double),
             cudaMemcpyDeviceToHost);
  CUERR

  status = 1;
  for (i = 0; i < test_grid.n_angular; ++i)
  {
    if (!(fabs(theta[i] - _at[i]) < eps))
    {
      printf("Fail! expected: %12.10f got: %12.10f\n", _at[i], theta[i]);
      status = 0;
    }
    if (!(fabs(phi[i] - _ap[i]) < eps))
    {
      printf("Fail! expected: %12.10f got: %12.10f\n", _ap[i], phi[i]);
      status = 0;
    }
    if (!(fabs(weights[i] - _aw[i]) < eps))
    {
      printf("Fail! expected: %12.10f got: %12.10f\n", _aw[i], weights[i]);
      status = 0;
    }
  }
  printf("Lebedev... %s\n", status ? "Passed" : "Failed");

  free(_at);
  free(_ap);
  free(_aw);
  free(theta);
  free(phi);
  free(weights);

  grid_free_cuda(&test_grid);
  grid_free(&grid);
}

int main(int argc, char const *argv[])
{
  test_becke_grid_init_cuda();
  return 0;
}