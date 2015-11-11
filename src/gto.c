/*file: gto.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "include/gto.h"

double gto_normalize_primitive(PrimitiveGaussian *f) {
  int i;
  double output;
  int aux = 0;

  for (i = 0; i < 3; ++i) {
    aux += f->l[i];
  }
  output = (pow((2.0 * f->exponent / M_PI), 0.75)) /
           sqrt(utils_factorial2(abs(2 * f->l[0] - 1)) *
                utils_factorial2(abs(2 * f->l[1] - 1)) *
                utils_factorial2(abs(2 * f->l[2] - 1)) /
                (pow((4.0 * f->exponent), aux)));

  return output;
}

void gto_compute_primitive(PrimitiveGaussian *f, double *coord, double *output,
                           const int n_coord) {
  int i, j, idx;
  double RP2 = 0.0, factor = 1.0, aux;

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, j, idx, aux, factor, RP2)
#endif
  for (i = 0; i < n_coord; ++i) {
    RP2 = 0.0;
    factor = 1.0;
    idx = i * 3;
    for (j = 0; j < 3; ++j) {
      aux = coord[idx + j] - f->origin[j];
      factor *= pow(aux, f->l[j]);
      RP2 += aux * aux;
    }
    output[i] =
        f->coefficient * f->normalization * factor * exp(-f->exponent * RP2);
  }
}

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

void gto_obaraSaika_recursion(double *x, double *y, double *z, double PA[3],
                              double PB[3], const double gamma, const int l_a,
                              const int l_b, const int max_l) {

  int i, j, idx, idxp1, idxm1;

  double pp;

  pp = 1.0 / (2.0 * gamma);

  x[0] = 1.0;
  y[0] = 1.0;
  z[0] = 1.0;

  x[1] = PB[0];
  y[1] = PB[1];
  z[1] = PB[2];

  i = 0;
  for (j = 1; j < l_b; j++) {
    idx = i * max_l + j;
    x[idx + 1] = PB[0] * x[j];
    y[idx + 1] = PB[1] * y[j];
    z[idx + 1] = PB[2] * z[j];
    x[idx + 1] += j * pp * x[idx - 1];
    y[idx + 1] += j * pp * y[idx - 1];
    z[idx + 1] += j * pp * z[idx - 1];
  }

  x[max_l] = PA[0];
  y[max_l] = PA[1];
  z[max_l] = PA[2];

  i = 1;
  for (j = 1; j <= l_b; j++) {
    idx = i * max_l + j;
    x[idx] = PA[0] * x[j];
    y[idx] = PA[1] * y[j];
    z[idx] = PA[2] * z[j];
    x[idx] += j * pp * x[j - 1];
    y[idx] += j * pp * y[j - 1];
    z[idx] += j * pp * z[j - 1];
  }

  for (i = 1; i < l_a; i++) {
    idx = i * max_l;
    idxp1 = (i + 1) * max_l;
    idxm1 = (i - 1) * max_l;
    x[idxp1] = PA[0] * x[idx];
    y[idxp1] = PA[1] * y[idx];
    z[idxp1] = PA[2] * z[idx];
    x[idxp1] += i * pp * x[idxm1];
    y[idxp1] += i * pp * y[idxm1];
    z[idxp1] += i * pp * z[idxm1];
    for (j = 1; j <= l_b; j++) {
      idx++;
      idxp1++;
      idxm1++;
      x[idxp1] = PA[0] * x[idx];
      y[idxp1] = PA[1] * y[idx];
      z[idxp1] = PA[2] * z[idx];
      x[idxp1] += i * pp * x[idxm1];
      y[idxp1] += i * pp * y[idxm1];
      z[idxp1] += i * pp * z[idxm1];
      x[idxp1] += j * pp * x[idx - 1];
      y[idxp1] += j * pp * y[idx - 1];
      z[idxp1] += j * pp * z[idx - 1];
    }
  }
}
