/*file: gauss_chebyshev.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#include "gauss_chebyshev.h"

void gaussChebyshev(int n, double* abscissas, double* weights)
{
  int i;
  double aux_a, aux_w, aux1, aux2;
  double aux = M_PI * 1.0 / (n + 1.0);

#ifdef _OMP
#pragma omp parallel for default(shared) private(i, aux_a, aux_w, aux1, aux2)
#endif
  for (i = 0; i < n; ++i)
  {
    aux_a = cos((i + 1) * aux);
    aux_w = (aux * (1.0 - aux_a * aux_a));

    /* Scale from interval (-1, 1) to (0, inf)*/
    /* This factor comes from the variable change r ! --> x */
    /* and from using chebyshev-gauss radial quadrature of second order */

    aux1 = (aux_a + 1.0) * 0.5;
    aux2 = aux1 * aux1;

    abscissas[i] = log(1.0 - (aux2 * aux2));
    weights[i] = ( (aux2 * aux1) / (sqrt(1.0 - aux_a * aux_a) * (1.0 - (aux2 * aux2) ) ) ) * aux_w * 2.0;
  }
}