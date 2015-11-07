/*file: system.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

extern "C" {
#include "include/system.h"
}

#include "include/cuda_helper.cuh"

void system_init_cuda(System *sys, System *sys_d) {
  int bytes_int = sys->n_particles * sizeof(int);
  int bytes_double = sys->n_particles * sizeof(double);

  /*Allocating space for the device structure*/
  cudaMalloc((void **)&sys_d->particle_number, bytes_int);
  cudaMalloc((void **)&sys_d->particle_radii, bytes_double);
  cudaMalloc((void **)&sys_d->particle_origin, 3 * bytes_double);

  /*Copying data to device*/
  sys_d->n_particles = sys->n_particles;
  cudaMemcpy(sys_d->particle_number, sys->particle_number, bytes_int,
             cudaMemcpyHostToDevice);
  cudaMemcpy(sys_d->particle_radii, sys->particle_radii, bytes_double,
             cudaMemcpyHostToDevice);
  cudaMemcpy(sys_d->particle_origin, sys->particle_origin, 3 * bytes_double,
             cudaMemcpyHostToDevice);

  /*Copying the BasisSet data*/
  basis_set_init_cuda(&sys->basis, &sys_d->basis);
}

void system_free_cuda(System *sys_d) {
  cudaFree(sys_d->particle_number);
  cudaFree(sys_d->particle_radii);
  cudaFree(sys_d->particle_origin);

  basis_set_free_cuda(&sys_d->basis);
}
