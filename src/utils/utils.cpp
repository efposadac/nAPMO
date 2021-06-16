/*
file: utils.cpp
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#include "utils.h"

int utils_factorial2(const int n) {

  if (n < 1) {
    return 0;
  } else {
    int output = 1;
    for (int i = n; i > 0; i -= 2) {
      output *= i;
    }
    return output;
  }
}

void utils_multiply_segmented_array(const int size, const int segments,
                                    double *array, double *output) {

  for (int i = 0; i < size; ++i) {
    double temp = array[i];
    for (int j = 1; j < segments; ++j) {
      temp *= array[j * size + i];
    }
    output[i] = temp;
  }
}

double utils_distance(double a[3], double b[3]) {

  double output = 0.0;
  for (int i = 0; i < 3; ++i) {
    double aux = a[i] - b[i];
    output += aux * aux;
  }

  return sqrt(output);
}