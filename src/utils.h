/*file: utils.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef UTILS_H
#define UTILS_H

#include <string.h>

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

static inline void multiply_segmented_array(int size, int segments, double* array, double * output){
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

#endif