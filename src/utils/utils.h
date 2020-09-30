/*file: utils.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

#ifndef UTILS_H
#define UTILS_H

#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

int utils_factorial2(const int n);

void utils_multiply_segmented_array(const int size, const int segments,
                                    double *array, double *output);

/*
Compute the Euclidian distance between two points.
*/
double utils_distance(double a[3], double b[3]);

#ifdef __cplusplus
}
#endif

#endif