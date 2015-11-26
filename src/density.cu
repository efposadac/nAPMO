/*file: density.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#define THREADS_PER_BLOCK 64

extern "C" {
#include "include/density.h"
}

void density_gto(BasisSet *basis, double *r, double *dens, double *output,
                 int size) {

  int sizeBasis;
  double *dens_d, *x_d, *y_d, *z_d, *output_d;
  BasisSet basis_d;

  /*
  Allocating Space
  */
  basis_set_init(basis, &basis_d);

  /*
  Copy points to device
  */
  cudaMalloc((void **)&x_d, size * sizeof(double));
  cudaMalloc((void **)&y_d, size * sizeof(double));
  cudaMalloc((void **)&z_d, size * sizeof(double));

  double *x = (double *)malloc(size * sizeof(double));
  double *y = (double *)malloc(size * sizeof(double));
  double *z = (double *)malloc(size * sizeof(double));

  for (int i = 0; i < size; ++i) {
    int idx = i * 3;
    x[i] = r[idx + 0];
    y[i] = r[idx + 1];
    z[i] = r[idx + 2];
  }

  cudaMemcpy(x_d, x, size * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(y_d, y, size * sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(z_d, z, size * sizeof(double), cudaMemcpyHostToDevice);

  free(x);
  free(y);
  free(z);

  /*
  Copy density matrix to device
  */
  sizeBasis = basis->n_cont * basis->n_cont;
  cudaMalloc((void **)&dens_d, sizeBasis * sizeof(double));
  cudaMemcpy(dens_d, dens, sizeBasis * sizeof(double), cudaMemcpyHostToDevice);

  /*
  Allocate space for output
  */
  cudaMalloc((void **)&output_d, size * sizeof(double));

  /*
  Calculate the function
  */
  dim3 dimGrid(((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1, 1);

  density_gto_kernel<<<dimGrid, THREADS_PER_BLOCK>>>(basis_d, x_d, y_d, z_d,
                                                     dens_d, output_d, size);
  CUERR

  /*
  Bring result back
  */
  cudaMemcpy(output, output_d, size * sizeof(double), cudaMemcpyDeviceToHost);

  /*
  clear memory
  */
  cudaFree(x_d);
  cudaFree(y_d);
  cudaFree(z_d);
  cudaFree(dens_d);
  cudaFree(output_d);
  basis_set_free(&basis_d);
}

__global__ void density_gto_kernel(BasisSet basis, double *x, double *y,
                                   double *z, double *dens, double *output,
                                   int size) {

  const unsigned int point = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  const unsigned int n_cont = basis.n_cont;

  double temp_val, function_value;
  double basis_val[64], r[3];

  if (point < size) {

    r[0] = x[point];
    r[1] = y[point];
    r[2] = z[point];

    basis_set_compute_gto(basis, r, basis_val);

    function_value = 0.0;
    for (int i = 0; i < n_cont; ++i) {

      temp_val = 0.0;
      for (int j = 0; j < n_cont; ++j) {
        temp_val += basis_val[j] * dens[i * n_cont + j];
      }

      function_value += temp_val * basis_val[i];
    }

    output[point] = function_value;
  }
}
