/*file: basis.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/
#ifndef BASIS_H
#define BASIS_H

struct _basis_set {
  int n_cont;            // Number of contractions.
  int *n_prim_cont;      // Number of primitives for each contraction.
  int *basis_l;          // Angular moment index sizeof (3*n_cont).
  double *exponent;      // exponent of each primitive.
  double *coefficient;   // contraction coefficient of each primitive.
  double *normalization; // Normalization constant for each primitive.
  double *origin;        // origin for each contraction.
};

typedef struct _basis_set BasisSet;

#ifdef _CUDA
/*
Copy the host BasisSet structure into the device.
*/
void basis_set_init_cuda(BasisSet *basis, BasisSet *basis_d);

/*
Free the memory used on the CUDA device.
*/
void basis_set_free_cuda(BasisSet *basis_d);

#endif

#endif
