/*file: system.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef SYSTEM_H
#define SYSTEM_H

#include "basis_set.h"

struct _system
{
  int n_particles;          // Number of particles in the system.
  int* particle_number;     // For atoms atomic number, else an identifier.
  double* particle_radii;   // Atomic radii, Particle radii?
  double* particle_origin;  // Origin of each atom / particle.
  BasisSet basis;
};
typedef struct _system System;

/*
CUDA functions
*/
#ifdef _CUDA

/*
Copy the host System structure into the device.
*/
void system_init_cuda(System* sys, System* sys_d);

/*
Free the memory used on the CUDA device.
*/
void system_free_cuda(System* sys_d);

#endif

#endif

/*double* particle_origin;  // Origin of each atom / particle.
  double* particle_radii;   // Atomic radii, Particle radii?
  int n_particles;          // Number of particles in the system.
  int* particle_number;     // For atoms atomic number, else an identifier.
  BasisSet basis;*/
