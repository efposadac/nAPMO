/*file: utils.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <string.h>

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

int utils_factorial2(const int n);

void utils_multiply_segmented_array(const int size, const int segments,
                                           double *array, double *output);

/*
Compute the Euclidian distance between two points.
*/
double utils_distance(double a[3], double b[3]);

#endif