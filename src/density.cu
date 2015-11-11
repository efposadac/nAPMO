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
  double *dens_d, *r_d, *output_d;
  BasisSet basis_d;

  /*
  Allocating Space
  */
  basis_set_init(basis, &basis_d);

  cudaMalloc((void **)&r_d, size * 3 * sizeof(double));
  cudaMemcpy(r_d, r, size * 3 * sizeof(double), cudaMemcpyHostToDevice);

  sizeBasis = basis->n_cont * basis->n_cont;
  cudaMalloc((void **)&dens_d, sizeBasis * sizeof(double));
  cudaMemcpy(dens_d, dens, sizeBasis * sizeof(double), cudaMemcpyHostToDevice);

  cudaMalloc((void **)&output_d, size * sizeof(double));

  /*
  Calculate the function
  */
  dim3 dimGrid(((size + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1, 1);

  density_gto_kernel<<<dimGrid, THREADS_PER_BLOCK>>>(basis_d, r_d, dens_d,
                                                     output_d, size);
  CUERR

  /*
  Bring result back
  */
  cudaMemcpy(output, output_d, size * sizeof(double), cudaMemcpyDeviceToHost);

  /*
  clear memory
  */
  cudaFree(r_d);
  cudaFree(dens_d);
  cudaFree(output_d);
  basis_set_free(&basis_d);
}

__global__ void density_gto_kernel(BasisSet basis, double *r, double *dens,
                                   double *output, int size) {

  const unsigned int point = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  const unsigned int n_cont = basis.n_cont;

  int i, j, idx, aux, counter;
  double temp_val, factor, RP2, temp, function_value;
  double basis_val[64];

  for (i = n_cont; i < n_cont * 2; ++i) {
    basis_val[i] = 0.0;
  }

  if (point < size) {
    idx = point * 3;
    counter = 0;
    temp_val = 0.0;

    for (j = 0; j < n_cont * 2; ++j) {
      basis_val[j] = 0.0;
    }

    for (i = 0; i < n_cont; ++i) {

      aux = i * 3;
      factor = 1.0, RP2 = 0.0;
      for (j = 0; j < 3; ++j) {
        temp = r[idx + j] - basis.origin[aux + j];
        RP2 += (temp * temp);
        factor *= pow(temp, basis.basis_l[aux + j]);
      }

      function_value = 0.0;
      for (j = 0; j < basis.n_prim_cont[i]; ++j) {
        function_value +=
            basis.coefficient[counter] * exp(-basis.exponent[counter] * RP2);
        counter += 1;
      }

      function_value *= factor * basis.normalization[i];

      for (j = 0; j < n_cont; ++j) {
        basis_val[n_cont + j] += function_value * dens[i * n_cont + j];
      }

      basis_val[i] = function_value;
    }

    for (i = 0; i < n_cont; ++i) {
      temp_val += basis_val[i] * basis_val[n_cont + i];
    }
    output[point] = temp_val;
  }
}
