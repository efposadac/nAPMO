/*file: utils.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "utils.h"

int utils_factorial2(const int n) {
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

void utils_multiply_segmented_array(const int size,
                                                  const int segments,
                                                  double *array,
                                                  double *output) {
  int i, j;
  double temp;

  for (i = 0; i < size; ++i) {
    temp = array[i];
    for (j = 1; j < segments; ++j) {
    temp *= array[j * size + i];
    }
    output[i] = temp;
  }
}

double utils_distance(double a[3], double b[3]) {
  int i;
  double output, aux;

  output = 0.0;
  for (i = 0; i < 3; ++i) {
    aux = a[i] - b[i];
    output += aux * aux;
  }

  return sqrt(output);
}