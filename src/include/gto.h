/*file: gto.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef GTO_H
#define GTO_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "utils.h"

struct _primitive {
  int l[3];
  double origin[3];
  double exponent;
  double coefficient;
  double normalization;
};
typedef struct _primitive PrimitiveGaussian;

/*
Calculates the normalization constant of this primitive.
*/
double gto_normalize_primitive(PrimitiveGaussian *f);

/*
Computes the value of the function at ``coord``.
*/
void gto_compute_primitive(PrimitiveGaussian *f, double *coord, double *output,
                           const int n_coord);

/*
Perform the Obara-Saika (1988) recursion to calculate the overlap integral:

:math:`<\phi_A|\phi_B>`

Where each :math:`\phi` corresponds to a GTO PrimitiveGaussian

Args:
    PA : Origin of :math:`\phi_A`
    PB : Origin of :math:`\phi_B`
    gamma : Reduced exponent (see: Gaussian product theorem)
    l_a : Angular moment index of :math:`\phi_A` + 2
    l_b : Angular moment index of :math:`\phi_B` + 2

Returns:
    x, y, z : x, y, and z components of the recursion.
*/
void gto_obaraSaika_recursion(double *x, double *y, double *z, double PA[3],
                              double PB[3], const double gamma, const int l_a,
                              const int l_b, const int max_l);

struct _contracted {
  int nprim;
  double normalization;
  double *origin;
  PrimitiveGaussian *prim;
};
typedef struct _primitive ContractedGaussian;

#endif