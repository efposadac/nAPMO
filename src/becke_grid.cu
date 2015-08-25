/*file: becke_grid_cuda.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "becke_grid.h"
}

#define THREADS_PER_BLOCK 8
#define THREADS_PER_BLOCK_2 64

/*
Implementation of functions, for documentation see the header file.
*/

void grid_init_cuda(Grid *grid, GridCuda *grid_d)
{
  grid_d->gridDim = make_int2(grid->n_radial, grid->n_angular);

  int bytes_radial = grid_d->gridDim.x * sizeof(double2);
  int bytes_angular = grid_d->gridDim.y * sizeof(double2);

  /* Allocate space for grid on device*/
  cudaMalloc((void **)&grid_d->rw, bytes_radial);
  cudaMalloc((void **)&grid_d->xy, bytes_angular);
  cudaMalloc((void **)&grid_d->zw, bytes_angular);

  /* Copying radial quadrature*/
  {
    double2 *buffer_rw = (double2 *)malloc(bytes_radial);

    for (int i = 0; i < grid_d->gridDim.x; ++i)
    {
      buffer_rw[i] = make_double2(grid->radial_abscissas[i], grid->radial_weights[i]);
    }

    cudaMemcpy(grid_d->rw, buffer_rw, bytes_radial, cudaMemcpyHostToDevice);

    free(buffer_rw);
  }

  /* Copying angular quadrature*/
  {
    double2 *buffer_xy = (double2 *)malloc(bytes_angular);
    double2 *buffer_zw = (double2 *)malloc(bytes_angular);

    for (int i = 0; i < grid_d->gridDim.y; ++i)
    {
      buffer_xy[i] = make_double2(grid->angular_theta[i], grid->angular_phi[i]);
      buffer_zw[i] = make_double2(1.0, grid->angular_weights[i]);
    }

    cudaMemcpy(grid_d->xy, buffer_xy, bytes_angular, cudaMemcpyHostToDevice);
    cudaMemcpy(grid_d->zw, buffer_zw, bytes_angular, cudaMemcpyHostToDevice);

    free(buffer_xy);
    free(buffer_zw);
  }
}

void grid_free_cuda(GridCuda *grid)
{
  cudaFree(grid->xy);
  cudaFree(grid->zw);
  cudaFree(grid->rw);
}

double grid_integrate_cuda(System *sys, Grid *grid)
{
  int sizeBasis, i;
  double *integral_d, integral;
  double *dens, *dens_d;

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

  for (i = 0; i < sizeBasis; ++i)
  {
    fscanf(file, "%lf", &dens[i]);
  }

  cudaMalloc((void **)&dens_d, sizeBasis * sizeof(double));
  cudaMemcpy(dens_d, dens, sizeBasis * sizeof(double), cudaMemcpyHostToDevice);
  free(dens);

  /*Initialize integral Value to 0.0*/
  cudaMalloc((void **)&integral_d, sizeof(double));
  cudaMemset(integral_d, 0.0, sizeof(double));

  /* Convert angular quadrature from spherical to cart*/
  {
    dim3 dimGrid(((grid_d.gridDim.y + THREADS_PER_BLOCK_2 - 1) / THREADS_PER_BLOCK_2), 1, 1);
    lebedev_to_cartesian_cuda <<<dimGrid, THREADS_PER_BLOCK_2>>>
        (grid_d.gridDim, grid_d.xy, grid_d.zw);
  }

  /* Perform integration*/
  {
    dim3 dimBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK, 1);
    dim3 dimGrid(((grid_d.gridDim.x + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK),
                 ((grid_d.gridDim.y + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1);

    grid_integrate_kernel <<<dimGrid, dimBlock>>>
        (sys_d, grid_d.gridDim, grid_d.xy, grid_d.zw, grid_d.rw, dens_d, integral_d);
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

__device__ double grid_weights_cuda(const int n_particles, double *__restrict__ particle_origin,
                                    double *__restrict__ particle_radii, double r[3], int particleID)
{
  double x_i[3], x_j, r_i, r_j, R_ij, mu_ij, rm_i, rm_j, chi, u_ij, a_ij, nu_ij;
  double sum = 0, output, aux1, aux2, aux3 = 0;
  int i, j, k, aux_i, aux_j;

  for (i = 0; i < n_particles; ++i)
  {
    aux1 = 1.0;
    aux_i = i * 3;
    x_i[0] = particle_origin[aux_i + 0];
    x_i[1] = particle_origin[aux_i + 1];
    x_i[2] = particle_origin[aux_i + 2];
    rm_i = particle_radii[i];

    for (j = 0; j < n_particles; ++j)
    {
      if (i != j)
      {
        // Internuclear distance (R_ij eq. 11)
        r_i = r_j = R_ij = 0.0;
        aux_j = j * 3;
        for (k = 0; k < 3; ++k)
        {
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
        if (fabs(a_ij) > 0.50)
        {
          a_ij *= 0.50 / fabs(a_ij);
        }
        // eq. A2
        nu_ij = mu_ij + a_ij * (1.0 - mu_ij * mu_ij);

        aux2 = grid_soft_mu_cuda(grid_soft_mu_cuda(grid_soft_mu_cuda(nu_ij)));
        aux1 *= 0.5 * (1.0 - aux2);
      }
    }
    sum += aux1;
    if (i == particleID) aux3 = aux1;
  }
  // eq. 22
  output = aux3 / sum;
  return output;
}

__device__ double grid_density_cuda(BasisSet basis, double *__restrict__ r, double *__restrict__ dens)
{
  int i, j, aux, counter = 0;
  int n_cont = basis.n_cont;
  double temp;
  double function_value, output = 0.0;
  double RP2, factor;
  double basis_val[64];

  for (i = n_cont; i < n_cont * 2; ++i)
  {
    basis_val[i] = 0.0;
  }

  for (i = 0; i < n_cont; ++i)
  {
    factor = 1.0, RP2 = 0.0;
    for (j = 0; j < 3; ++j)
    {
      aux = i * 3;
      temp = r[j] - basis.origin[aux + j];
      RP2 += (temp * temp);
      factor *= pow(temp, basis.basis_l[aux + j]);
    }

    function_value = 0.0;
    for (j = 0; j < basis.n_prim_cont[i]; ++j)
    {
      function_value += basis.coefficient[counter] * exp(-basis.exponent[counter] * RP2);
      counter += 1;
    }

    function_value *= factor * basis.normalization[i];

    for (j = 0; j < n_cont; ++j)
    {
      basis_val[n_cont + j] += function_value * dens[i * n_cont + j];
    }

    basis_val[i] = function_value;
  }

  for (i = 0; i < n_cont; ++i)
  {
    output += basis_val[i] * basis_val[n_cont + i];
  }

  return output;
}

__global__ void grid_integrate_kernel(const System sys, const int2 gridDim, double2 *__restrict__ xy,
                                      double2 *__restrict__ zw, double2 *__restrict__ rw,
                                      double *__restrict__ dens, double *__restrict__ integral)
{
  short k, aux_k;
  double rm, rad, p, f, r[3];

  __shared__ double temp[THREADS_PER_BLOCK_2];
  __shared__ double sum_block;

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  const unsigned int j = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;
  const unsigned int l = __umul24(threadIdx.x, THREADS_PER_BLOCK) + threadIdx.y;

  temp[l] = 0.0;
  sum_block = 0.0;

  if (i < gridDim.x && j < gridDim.y)
  {
    /* Calculate r*/
    const double2 aux_rw = rw[i];
    const double2 aux_xy = xy[j];
    const double2 aux_zw = zw[j];

    const double auxFactor = aux_rw.y * aux_zw.y;

    for (k = 0; k < sys.n_particles; ++k)
    {
      rm = sys.particle_radii[k];

      if (sys.particle_number[k] != 1)
      {
        rm *= 0.5;
      }

      rad = rm * aux_rw.x;

      aux_k = k * 3;
      r[0] = rad * aux_xy.x + sys.particle_origin[aux_k + 0];
      r[1] = rad * aux_xy.y + sys.particle_origin[aux_k + 1];
      r[2] = rad * aux_zw.x + sys.particle_origin[aux_k + 2];

      /*Calculate Becke weights */
      p = grid_weights_cuda(sys.n_particles, sys.particle_origin, sys.particle_radii, r, k);

      /*Calculate the functional*/
      f = grid_density_cuda(sys.basis, r, dens);

      /*Calculate Integral */
      temp[l] += (rad * rad * rm * auxFactor * p * f);
    }
    __syncthreads();
    atomicAdd(&sum_block, temp[l]);
  }
  __syncthreads();
  if (l == 0)
  {
    atomicAdd(integral, sum_block);
  }
}
