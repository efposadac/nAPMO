/*file: becke_grid_cuda.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "include/becke_grid.h"
}

#define THREADS_PER_BLOCK 64

/*
Implementation of functions, for documentation see the header file.
*/

void grid_init_cuda(BeckeGrid *grid, GridCuda *grid_d) {
  grid_d->gridDim = make_int2(grid->size, grid->ncenter);

  /* Allocate space for grid on device*/
  int ncenter = grid_d->gridDim.y * sizeof(double2);
  cudaMalloc((void **)&grid_d->radii, ncenter);
  cudaMalloc((void **)&grid_d->origin, ncenter * 3);

  int size = grid_d->gridDim.x * sizeof(double2);
  cudaMalloc((void **)&grid_d->xy, size);
  cudaMalloc((void **)&grid_d->zw, size);

  /* Copying atomic information*/
  {
    cudaMemcpy(grid_d->radii, grid->radii, ncenter, cudaMemcpyHostToDevice);
    cudaMemcpy(grid_d->origin, grid->origin, ncenter * 3,
               cudaMemcpyHostToDevice);
  }

  /* Copying grid points*/
  {
    double2 *buffer_xy = (double2 *)malloc(size);
    double2 *buffer_zw = (double2 *)malloc(size);

    for (int i = 0; i < grid_d->gridDim.x; ++i) {
      int idx = i * 3;
      buffer_xy[i] = make_double2(grid->points[idx], grid->points[idx + 1]);
      buffer_zw[i] = make_double2(grid->points[idx + 2], grid->weights[i]);
    }

    cudaMemcpy(grid_d->xy, buffer_xy, size, cudaMemcpyHostToDevice);
    cudaMemcpy(grid_d->zw, buffer_zw, size, cudaMemcpyHostToDevice);

    free(buffer_xy);
    free(buffer_zw);
  }
}

void grid_free_cuda(GridCuda *grid) {
  cudaFree(grid->radii);
  cudaFree(grid->origin);
  cudaFree(grid->xy);
  cudaFree(grid->zw);
}

double grid_integrate_cuda(System *sys, BeckeGrid *grid, double *rad,
                           int nrad) {
  int sizeBasis, i;
  double *integral_d, integral;
  double *dens, *dens_d, *rad_d;

  /*Initialize system data in the device*/
  System sys_d;
  system_init_cuda(sys, &sys_d);

  /* Initialize grid on device*/
  GridCuda grid_d;
  grid_init_cuda(grid, &grid_d);

  /*Fetch density file and copy to device memory.*/
  sizeBasis = sys->basis.n_cont * sys->basis.n_cont;
  dens = (double *)malloc(sizeBasis * sizeof(double));

  /* TODO: this file has to disappear*/
  FILE *file = fopen("data.dens", "r");

  for (i = 0; i < sizeBasis; ++i) {
    fscanf(file, "%lf", &dens[i]);
  }

  cudaMalloc((void **)&dens_d, sizeBasis * sizeof(double));
  cudaMemcpy(dens_d, dens, sizeBasis * sizeof(double), cudaMemcpyHostToDevice);
  free(dens);

  cudaMalloc((void **)&rad_d, nrad * grid_d.gridDim.y * sizeof(double));
  cudaMemcpy(rad_d, rad, nrad * grid_d.gridDim.y * sizeof(double), cudaMemcpyHostToDevice);

  /*Initialize integral Value to 0.0*/
  cudaMalloc((void **)&integral_d, sizeof(double));
  cudaMemset(integral_d, 0.0, sizeof(double));

  /* Perform integration*/
  {
    dim3 dimGrid(
        (((grid_d.gridDim.x/grid_d.gridDim.y) + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1,
        1);

    grid_integrate_kernel<<<dimGrid, THREADS_PER_BLOCK>>>(
        sys_d, grid_d.gridDim, grid_d.radii, grid_d.origin, grid_d.xy,
        grid_d.zw, dens_d, rad_d, nrad, integral_d);
    CUERR
  }

  /* Bring result to device */
  integral = 0.0;
  cudaMemcpy(&integral, integral_d, sizeof(double), cudaMemcpyDeviceToHost);

  /* Free memory */
  cudaFree(dens_d);
  cudaFree(integral_d);

  system_free_cuda(&sys_d);
  grid_free_cuda(&grid_d);

  return integral * 4.0 * M_PI;
}

/*
CUDA kernels:
*/

__global__ void
grid_integrate_kernel(const System sys, const int2 gridDim,
                      double *__restrict__ radii, double *__restrict__ origin,
                      double2 *__restrict__ xy, double2 *__restrict__ zw,
                      double *__restrict__ dens, double *__restrict__ rad,
                      int nrad, double *__restrict__ integral) {
  short k;
  double p, f, r[3];

  __shared__ double temp[THREADS_PER_BLOCK];
  __shared__ double sum_block;

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  temp[threadIdx.x] = 0.0;
  sum_block = 0.0;
  int size = gridDim.x / gridDim.y;
  int size2 = size / nrad;

  if (i < size) {
    for (k = 0; k < gridDim.y; ++k) {
    /* Calculate r*/
    const double2 aux_xy = xy[k * size + i];
    const double2 aux_zw = zw[k * size + i];

      /*Calculate Becke weights */
      r[0] = aux_xy.x;
      r[1] = aux_xy.y;
      r[2] = aux_zw.x;

      p = grid_weights_cuda(gridDim.y, origin, radii, r, k);

      /*Calculate the functional*/
      f = grid_density_cuda(sys.basis, r, dens);

      /*Calculate Integral */
      int idxr = k * nrad + __float2int_rz(i / size2);
      temp[threadIdx.x] += (rad[idxr] * rad[idxr] * radii[k] * aux_zw.y * p * f);
    }
    __syncthreads();
    atomicAdd(&sum_block, temp[threadIdx.x]);
  }
  __syncthreads();
  if (threadIdx.x == 0) {
    atomicAdd(integral, sum_block);
  }
}

__device__ double grid_weights_cuda(const int n_particles,
                                    double *__restrict__ particle_origin,
                                    double *__restrict__ particle_radii,
                                    double r[3], int particleID) {
  double x_i[3], x_j, r_i, r_j, R_ij, mu_ij, rm_i, rm_j, chi, u_ij, a_ij, nu_ij;
  double sum = 0, output, aux1, aux2, aux3 = 0;
  int i, j, k, aux_i, aux_j;

  for (i = 0; i < n_particles; ++i) {
    aux1 = 1.0;
    aux_i = i * 3;
    x_i[0] = particle_origin[aux_i + 0];
    x_i[1] = particle_origin[aux_i + 1];
    x_i[2] = particle_origin[aux_i + 2];
    rm_i = particle_radii[i];

    for (j = 0; j < n_particles; ++j) {
      if (i != j) {
        // Internuclear distance (R_ij eq. 11)
        r_i = r_j = R_ij = 0.0;
        aux_j = j * 3;
        for (k = 0; k < 3; ++k) {
          x_j = particle_origin[aux_j + k];
          r_i += (r[k] - x_i[k]) * (r[k] - x_i[k]);
          r_j += (r[k] - x_j) * (r[k] - x_j);
          R_ij += (x_i[k] - x_j) * (x_i[k] - x_j);
        }

        r_i = sqrt(r_i);
        r_j = sqrt(r_j);
        R_ij = sqrt(R_ij);

        // \mu_ij eq. 11
        mu_ij = (r_i - r_j) / R_ij;

        // Atomic size adjustment. see appendix, Becke, 1988.
        rm_j = particle_radii[j];

        // eq. A4
        chi = rm_i / rm_j;
        // eq. A6
        u_ij = (chi - 1.0) / (chi + 1.0);
        // eq. A5
        a_ij = u_ij / ((u_ij * u_ij) - 1.0);
        // eq. A3
        if (fabs(a_ij) > 0.50) {
          a_ij *= 0.50 / fabs(a_ij);
        }
        // eq. A2
        nu_ij = mu_ij + a_ij * (1.0 - mu_ij * mu_ij);

        aux2 = grid_soft_mu_cuda(grid_soft_mu_cuda(grid_soft_mu_cuda(nu_ij)));
        aux1 *= 0.5 * (1.0 - aux2);
      }
    }
    sum += aux1;
    if (i == particleID)
      aux3 = aux1;
  }
  // eq. 22
  output = aux3 / sum;
  return output;
}

__device__ double grid_density_cuda(BasisSet basis, double *__restrict__ r,
                                    double *__restrict__ dens) {
  int i, j, aux, counter = 0;
  int n_cont = basis.n_cont;
  double temp;
  double function_value, output = 0.0;
  double RP2, factor;
  double basis_val[64];

  for (i = n_cont; i < n_cont * 2; ++i) {
    basis_val[i] = 0.0;
  }

  for (i = 0; i < n_cont; ++i) {
    factor = 1.0, RP2 = 0.0;
    for (j = 0; j < 3; ++j) {
      aux = i * 3;
      temp = r[j] - basis.origin[aux + j];
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
    output += basis_val[i] * basis_val[n_cont + i];
  }

  return output;
}
