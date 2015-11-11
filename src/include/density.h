/*file: density.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef DENSITY_H
#define DENSITY_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "basis_set.h"

#ifdef _OMP
#include <omp.h>
#endif

/*
Calculate the density at point r
*/
void density_gto(BasisSet *basis, double *r, double *dens, double *output,
                 int size);

#ifdef _CUDA
#include "cuda_helper.cuh"

/*
Calculate the density at point r (CUDA version)
*/
__global__ void density_gto_kernel(BasisSet basis, double *r, double *dens,
                                   double *output, int size);
#endif

#endif