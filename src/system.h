/*file: system.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef SYSTEM_H
#define SYSTEM_H

#include "basis_set.h"

struct _system {
  int n_particles;          // Number of particles in the system.
  int* particle_number;     // For atoms atomic number, else an identifier.
  double* particle_radii;   // Atomic radii, Particle radii?
  double* particle_origin;  // Origin of each atom / particle.
  BasisSet basis;
};
typedef struct _system System;

#endif
