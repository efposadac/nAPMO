/*file: utils.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <math.h>

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

static inline int utils_factorial2(const int n) {
  int i;
  int output = 1;

  if (n < 1) {
    return 0;
  } else {
    for (i = n; i > 0; i -= 2) {
      output *= i;
    }
    return output;
  }
}

static inline void utils_multiply_segmented_array(const int size,
                                                  const int segments,
                                                  double *array,
                                                  double *output) {
  int i, j;
  memcpy(output, array, size * sizeof(double));

  for (i = 1; i < segments; ++i) {
#ifdef _OMP
#pragma omp parallel for default(shared) firstprivate(i) private(j)
#endif
    for (j = 0; j < size; ++j) {

      output[j] *= array[i * size + j];
    }
  }
}

/*
Compute the Euclidian distance between two points.
*/
static inline double utils_distance(double a[3], double b[3]) {
  int i;
  double output, aux;

  output = 0.0;
  for (i = 0; i < 3; ++i) {
    aux = a[i] - b[i];
    output += aux * aux;
  }

  return sqrt(output);
}

#endif