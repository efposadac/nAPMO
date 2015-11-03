/*file: grid.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef GRID_H
#define GRID_H

#include "lebedev.h"

void grid_evaluate_atomic_expansion(int lmax, int lorder, int size,
                                        double *decomposition, double *points,
                                        double *output);

double grid_atomic_integrate(const int size, const int segments, double rm, double *f,
                         double *w);
#endif