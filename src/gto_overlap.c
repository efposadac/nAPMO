/*file: gto_overlap.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "include/gto.h"


double gto_overlap_primitive(PrimitiveGaussian *f_a, PrimitiveGaussian *f_b) {
  int i;
  int l_a = 0, l_b = 0, max_l;

  double output;
  double AB, P0;
  double PA[3], PB[3], AB2 = 0.0;
  double gamma, gammaInv, preFactor;

  double x0, y0, z0;
  double *x, *y, *z;

  gamma = f_a->exponent + f_b->exponent;
  gammaInv = 1.0 / gamma;

  for (i = 0; i < 3; ++i) {
    AB = f_a->origin[i] - f_b->origin[i];
    P0 = (f_a->exponent * f_a->origin[i] + f_b->exponent * f_b->origin[i]) *
         gammaInv;
    PA[i] = P0 - f_a->origin[i];
    PB[i] = P0 - f_b->origin[i];
    AB2 += AB * AB;
    l_a += f_a->l[i];
    l_b += f_b->l[i];
  }

  preFactor = exp(-f_a->exponent * f_b->exponent * AB2 * gammaInv) *
              (sqrt(M_PI * gammaInv) * M_PI * gammaInv * f_a->coefficient *
               f_b->coefficient * f_a->normalization * f_b->normalization);

  // recursion
  max_l = max(l_a, l_b) + 4;

  x = (double *)calloc(max_l * max_l, sizeof(double));
  y = (double *)calloc(max_l * max_l, sizeof(double));
  z = (double *)calloc(max_l * max_l, sizeof(double));

  gto_obaraSaika_recursion(x, y, z, PA, PB, gamma, l_a + 2, l_b + 2, max_l);

  // Calculating integrals for primitives
  x0 = x[f_a->l[0] * max_l + f_b->l[0]];
  y0 = y[f_a->l[1] * max_l + f_b->l[1]];
  z0 = z[f_a->l[2] * max_l + f_b->l[2]];

  output = preFactor * x0 * y0 * z0;

  free(x);
  free(y);
  free(z);

  return output;
}