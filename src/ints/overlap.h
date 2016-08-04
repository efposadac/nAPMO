/*file: overlap.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef OVERLAP_H
#define OVERLAP_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../system/gto.h"

#ifdef __cplusplus
extern "C" {
#endif

class PrimitiveGaussian;

/*
Computes the overlap integral between two primitive Gaussian functions.
Used only for normalization!
*/
double overlap_primitive(PrimitiveGaussian *f_a, PrimitiveGaussian *f_b);

#ifdef __cplusplus
}
#endif

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
void obaraSaika_recursion(double *x, double *y, double *z, double PA[3],
                          double PB[3], const double gamma, const int l_a,
                          const int l_b, const int max_l);

#endif