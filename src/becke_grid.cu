/*file: becke_grid_cuda.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "becke_grid.h"
#include "cuda_helper.cuh"
}

#define THREADS_PER_BLOCK 8

void grid_init_cuda(Grid *grid, Grid *grid_d)
{
  double bytes_radial = grid->n_radial * sizeof(double);
  double bytes_angular = grid->n_angular * sizeof(double);

  /*Allocating space for the device structure*/
  cudaMalloc((void **)&grid_d->radial_abscissas, bytes_radial);
  cudaMalloc((void **)&grid_d->radial_weights, bytes_radial);

  cudaMalloc((void **)&grid_d->angular_theta, bytes_angular);
  cudaMalloc((void **)&grid_d->angular_phi, bytes_angular);
  cudaMalloc((void **)&grid_d->angular_weights, bytes_angular);

  /*Copying data to device*/
  grid_d->n_radial = grid->n_radial;
  grid_d->n_angular = grid->n_angular;
  cudaMemcpy(grid_d->angular_theta, grid->angular_theta, bytes_angular, cudaMemcpyHostToDevice);
  cudaMemcpy(grid_d->angular_phi, grid->angular_phi, bytes_angular, cudaMemcpyHostToDevice);
  cudaMemcpy(grid_d->angular_weights, grid->angular_weights, bytes_angular, cudaMemcpyHostToDevice);

  /* Calculate Gauss-Chebyshev points in device (not copying)*/
  gaussChebyshev_cuda(grid->n_radial, grid_d->radial_abscissas, grid_d->radial_weights);
}

void grid_free_cuda(Grid *grid)
{
  cudaFree(grid->radial_abscissas);
  cudaFree(grid->radial_weights);
  cudaFree(grid->angular_theta);
  cudaFree(grid->angular_phi);
  cudaFree(grid->angular_weights);
}

double grid_integrate_cuda(System *sys, Grid *grid)
{
  int i;
  double *integral_d, integral;

  /*Initialize system data in the device*/
  System sys_d;
  system_init_cuda(sys, &sys_d);

  /* Initialize grid on device*/
  Grid grid_d;
  grid_init_cuda(grid, &grid_d);

  /*Fetch density file and copy to device memory.*/
  int size = sys->basis.n_cont * sys->basis.n_cont;
  double *dens = (double *)malloc(size * sizeof(double));

  FILE *file;
  file = fopen("data.dens", "r");

  for (i = 0; i < size; ++i)
  {
    fscanf(file, "%lf", &dens[i]);
  }

  double *dens_d;
  cudaMalloc((void **)&dens_d, size * sizeof(double));
  cudaMemcpy(dens_d, dens, size * sizeof(double), cudaMemcpyHostToDevice);
  free(dens);

  /* Initialize integral Value to 0*/
  cudaMalloc((void **)&integral_d, sizeof(double));
  cudaMemset(integral_d, 0, sizeof(double));

  /*Run integrator*/
  dim3 dimBlock(THREADS_PER_BLOCK, THREADS_PER_BLOCK, 1);

  dim3 dimGrid(((grid_d.n_radial + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK),
               ((grid_d.n_angular + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK), 1);

  grid_integrate_kernel <<<dimGrid, dimBlock>>> (sys_d, grid_d, dens_d, integral_d);
  CUERR

  /* Bring result to device */
  integral = 0.0;
  cudaMemcpy(&integral, integral_d, sizeof(double), cudaMemcpyDeviceToHost);

  // /* Free memory */
  cudaFree(dens_d);
  cudaFree(integral_d);

  system_free_cuda(&sys_d);
  grid_free_cuda(&grid_d);

  return integral * 8.0 * M_PI;
}

__device__ double atomicAdd(double *address, double val)
{
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do
  {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

__device__ double grid_soft_mu_cuda(const double mu) { return 0.5 * mu * (3.0 - mu * mu); }

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

__device__ double grid_density_cuda(BasisSet basis, double * __restrict__ r, double * __restrict__ dens)
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

__global__ void grid_integrate_kernel(System sys, Grid grid, double * __restrict__ dens, double *integral)
{
  double rm, rad, r[3], factor;
  double aux1, aux2, aux3, aux4, auxLog, auxFactor, auxSqrt, aux1m4, auxSqrtx1m4;
  double sin_t, cos_t, sin_p, cos_p;
  double q_r, w_r, w_a, p, f;
  double sum;
  int k, auxK;

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;
  const unsigned int j = __umul24(blockIdx.y, blockDim.y) + threadIdx.y;

  if ((i < grid.n_radial) && (j < grid.n_angular))
  {
    q_r = grid.radial_abscissas[i];
    w_r = grid.radial_weights[i];
    w_a = grid.angular_weights[j];

    aux1 = q_r * 0.5 + 0.5;
    aux2 = aux1 * aux1;
    aux3 = aux2 * aux1;
    aux4 = aux2 * aux2;
    auxLog = log(1.0 - aux4);
    auxFactor = w_r * w_a * aux3;
    auxSqrt = sqrt(1.0 - q_r * q_r);
    aux1m4 = 1.0 - aux4;
    auxSqrtx1m4 = auxSqrt * aux1m4;
    auxFactor = auxFactor / auxSqrtx1m4;

    sum = 0.0;

    sincos(grid.angular_theta[j], &sin_t, &cos_t);
    sincos(grid.angular_phi[j], &sin_p, &cos_p);

    for (k = 0; k < sys.n_particles; ++k)
    {
      rm = sys.particle_radii[k];

      if (sys.particle_number[k] != 1)
      {
        rm *= 0.5;
      }

      rad = -rm * auxLog;

      factor = rad * rad * rm * auxFactor;

      auxK = k * 3;
      r[0] = (rad * sin_t * cos_p) + sys.particle_origin[auxK + 0];
      r[1] = (rad * sin_t * sin_p) + sys.particle_origin[auxK + 1];
      r[2] = (rad * cos_t) + sys.particle_origin[auxK + 2];

      /*Calculate Becke weights */
      p = grid_weights_cuda(sys.n_particles, sys.particle_origin, sys.particle_radii, r, k);

      /*Calculate the functional*/
      f = grid_density_cuda(sys.basis, r, dens);

      /*Calculate Integral */
      sum += (factor * p * f);
    }
    atomicAdd(integral, sum);
  }
}
