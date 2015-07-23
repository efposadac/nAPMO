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

#define THREADS_PER_BLOCK 64

/*
Iterated cutoff profile. eq. 21, Becke 1988. (CUDA Device version)

Note:
    Avoid the use of __host__ __device__ in order to allow the compilation with
other compilers
    in the case of non CUDA compilation.

*/
__device__ __forceinline__ double grid_soft_mu_cuda(const double &mu)
{
  return 0.5 * mu * (3.0 - mu * mu);
}

/*
Implementation of double2 multiply operation.
*/
__device__ __forceinline__ double2 mult_double2(const double2 &a, const double2 &b)
{
  double2 r;
  r.x = a.x * b.x;
  r.y = a.y * b.y;
  return r;
}

/*
Implementation of atomicAdd for double precision variables.
*/
__device__ __forceinline__ double atomicAdd(double *address, double val)
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

/*
Implementation of functions, for documentation see the header file.
*/

void grid_init_cuda(Grid *grid, GridCuda *grid_d)
{
  int bytes_grid, dimGrid;
  ;
  double *data;
  double2 *data2, *buffer2;

  grid_d->gridDim = make_int2(grid->n_radial, grid->n_angular);
  bytes_grid = grid_d->gridDim.x * grid_d->gridDim.y * sizeof(double2);

  /* Allocate space for grid on device*/
  cudaMalloc((void **)&grid_d->xy, bytes_grid);
  cudaMalloc((void **)&grid_d->zr, bytes_grid);
  cudaMalloc((void **)&grid_d->wf, bytes_grid);

  /* Convert grid coordinates to cartesian coordinates */

  /* Angular spherical coordinates to cartesian with r = 1*/
  cudaMalloc((void **)&data2, grid_d->gridDim.y * sizeof(double2));
  cudaMalloc((void **)&data, grid_d->gridDim.y * sizeof(double));

  buffer2 = (double2 *)malloc(grid_d->gridDim.y * sizeof(double2));

  for (int i = 0; i < grid_d->gridDim.y; ++i)
  {
    buffer2[i] = make_double2(grid->angular_theta[i], grid->angular_phi[i]);
  }

  cudaMemcpy(data2, buffer2, grid_d->gridDim.y * sizeof(double2), cudaMemcpyHostToDevice);
  cudaMemcpy(data, grid->angular_weights, grid_d->gridDim.y * sizeof(double), cudaMemcpyHostToDevice);

  dimGrid = ((grid_d->gridDim.y + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  grid_ang_sph_to_cart_kernel << <dimGrid, THREADS_PER_BLOCK>>>
      (grid_d->gridDim, grid_d->xy, grid_d->zr, grid_d->wf, data2, data);
  CUERR

  free(buffer2);
  cudaFree(data2);
  cudaFree(data);

  /* Complete conversion calculating r */
  cudaMalloc((void **)&data2, grid_d->gridDim.x * sizeof(double2));

  buffer2 = (double2 *)malloc(grid_d->gridDim.x * sizeof(double2));

  for (int i = 0; i < grid_d->gridDim.x; ++i)
  {
    buffer2[i] = make_double2(grid->radial_abscissas[i], grid->radial_weights[i]);
  }

  cudaMemcpy(data2, buffer2, grid_d->gridDim.x * sizeof(double2), cudaMemcpyHostToDevice);

  dimGrid = ((grid_d->gridDim.y + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);
  grid_rad_sph_to_cart_kernel << <dimGrid, THREADS_PER_BLOCK>>>
      (grid_d->gridDim, grid_d->xy, grid_d->zr, grid_d->wf, data2);
  CUERR

  free(buffer2);
  cudaFree(data2);
}

void grid_free_cuda(GridCuda *grid)
{
  cudaFree(grid->xy);
  cudaFree(grid->zr);
  cudaFree(grid->wf);
}

double grid_integrate_cuda(System *sys, Grid *grid)
{
  int dimGrid, sizeGrid, sizeBasis, i;
  double *integral_d, integral;
  double *dens, *dens_d;

  FILE *file;

  /*Initialize system data in the device*/
  System sys_d;
  system_init_cuda(sys, &sys_d);

  /* Initialize grid on device*/
  GridCuda grid_d;
  grid_init_cuda(grid, &grid_d);

  /*Fetch density file and copy to device memory.*/
  sizeBasis = sys->basis.n_cont * sys->basis.n_cont;
  dens = (double *)malloc(sizeBasis * sizeof(double));

  file = fopen("data.dens", "r");

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

  /* Perform integration*/
  sizeGrid = grid_d.gridDim.x * grid_d.gridDim.y;
  dimGrid = ((sizeGrid + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK);

  grid_integrate_kernel <<<dimGrid, THREADS_PER_BLOCK>>>
      (sys_d, grid_d.gridDim, grid_d.xy, grid_d.zr, grid_d.wf, dens_d, integral_d);
  CUERR

  /* Bring result to device */
  integral = 0.0;
  cudaMemcpy(&integral, integral_d, sizeof(double), cudaMemcpyDeviceToHost);

  /* Free memory */
  cudaFree(dens_d);
  cudaFree(integral_d);

  system_free_cuda(&sys_d);
  grid_free_cuda(&grid_d);

  return integral * 8.0 * M_PI;
}

/*
CUDA kernels:
*/

__global__ void grid_ang_sph_to_cart_kernel(const int2 gridDim, double2 *__restrict__ xy,
                                            double2 *__restrict__ zr, double2 *__restrict__ wf,
                                            double2 *__restrict__ data_tp, double *__restrict__ data_w)
{
  __shared__ double2 aux_wf[THREADS_PER_BLOCK];
  __shared__ double2 aux_xy[THREADS_PER_BLOCK];
  __shared__ double2 aux_zr[THREADS_PER_BLOCK];

  unsigned int ij;
  double sin_t, sin_p, cos_t, cos_p;

  const unsigned int j = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  if (j < gridDim.y)
  {
    const double2 tp = data_tp[j];

    sincos(tp.x, &sin_t, &cos_t);
    sincos(tp.y, &sin_p, &cos_p);

    const double aux1 = sin_t * cos_p;
    const double aux2 = sin_t * sin_p;

    aux_xy[threadIdx.x] = make_double2(aux1, aux2);
    aux_zr[threadIdx.x] = make_double2(cos_t, 1.0);
    aux_wf[threadIdx.x] = make_double2(data_w[j], 1.0);

#pragma unroll
    for (int i = 0; i < gridDim.x; ++i)
    {
      ij = __umul24(i, gridDim.y) + j;
      xy[ij] = aux_xy[threadIdx.x];
      zr[ij] = aux_zr[threadIdx.x];
      wf[ij] = aux_wf[threadIdx.x];
    }
  }
}

__global__ void grid_rad_sph_to_cart_kernel(const int2 gridDim, double2 *__restrict__ xy,
                                            double2 *__restrict__ zr, double2 *__restrict__ wf,
                                            double2 *__restrict__ data_rw)
{
  __shared__ double2 aux_log[THREADS_PER_BLOCK];
  __shared__ double2 aux_wf[THREADS_PER_BLOCK];

  unsigned int ij;
  double aux1, aux2, aux3, aux4, auxLog;
  double auxSqrt, auxSqrtx1m4, auxFactor;

  const unsigned int j = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  if (j < gridDim.y)
  {
#pragma unroll
    for (int i = 0; i < gridDim.x; ++i)
    {
      const double2 rw = data_rw[i];
      aux1 = rw.x * 0.5 + 0.5;
      aux2 = aux1 * aux1;
      aux3 = aux1 * aux2;
      aux4 = aux2 * aux2;
      auxLog = log(1.0 - aux4);

      auxSqrt = __dsqrt_rd(__fma_rd(-rw.x, rw.x, 1.0));
      auxSqrtx1m4 = __fma_rd(-auxSqrt, aux4, auxSqrt);
      auxFactor = __drcp_rd(auxSqrtx1m4);

      aux_log[threadIdx.x] = make_double2(auxLog, auxLog);
      aux_wf[threadIdx.x] = make_double2(rw.y * aux3, auxFactor);

      ij = __umul24(i, gridDim.y) + j;
      xy[ij] = mult_double2(xy[ij], aux_log[threadIdx.x]);
      zr[ij] = mult_double2(zr[ij], aux_log[threadIdx.x]);
      wf[ij] = mult_double2(wf[ij], aux_wf[threadIdx.x]);
    }
  }
}

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
                                      double2 *__restrict__ zr, double2 *__restrict__ wf,
                                      double *__restrict__ dens, double *__restrict__ integral)
{
  short k, aux_k;
  double rm, rad, factor, sum, p, f, r[3];
  double2 aux_xy, aux_zr, aux_wf;

  const unsigned int i = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

  sum = 0.0;
  if (i < __umul24(gridDim.x, gridDim.y))
  {
    aux_xy = xy[i];
    aux_zr = zr[i];
    aux_wf = wf[i];

    for (k = 0; k < sys.n_particles; ++k)
    {
      rm = sys.particle_radii[k];

      if (sys.particle_number[k] != 1)
      {
        rm *= 0.5;
      }

      rad = -rm * aux_zr.y;
      factor = rad * rad * rm * aux_wf.y * aux_wf.x;

      aux_k = k * 3;
      r[0] = -rm * aux_xy.x + sys.particle_origin[aux_k + 0];
      r[1] = -rm * aux_xy.y + sys.particle_origin[aux_k + 1];
      r[2] = -rm * aux_zr.x + sys.particle_origin[aux_k + 2];

      /*Calculate Becke weights */
      p = grid_weights_cuda(sys.n_particles, sys.particle_origin, sys.particle_radii, r, k);

      /*Calculate the functional*/
      f = grid_density_cuda(sys.basis, r, dens);

      /*Calculate Integral */
      sum += (factor * p * f);
    }
    // *integral = 1.0;
    atomicAdd(integral, sum);
  }
}
