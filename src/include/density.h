/*file: density.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

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
Calculate the density values for an array of points r
*/
void density_gto(BasisSet *basis, double *r, double *dens, double *output,
                 int size);

/*
Calculate the density at point r
*/
double density_gto_r(BasisSet *basis, double r[3], double *dens);

#ifdef _CUDA
#include "cuda_helper.cuh"

/*
Calculate the density at point r (CUDA version)
*/
__global__ void density_gto_kernel(BasisSet basis, double *x, double *y,
                                   double *z, double *dens, double *output,
                                   double *basis_val, int size);
#endif

#endif