/*file: radial.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef RADIAL_H
#define RADIAL_H

#include "../horton/rtransform.h"
#include "../utils/utils.h"
#include <math.h>
#include <stdlib.h>

struct RadialGrid {

private:
  unsigned int size; // number of grid points.
  double radii;      // covalent radius of atom
  double *points;    // points of the grid
  double *weights;   // weights of the grid
  RTransform *rtf;

public:
  RadialGrid() : size(0), radii(1.0){};

  RadialGrid(const RadialGrid &) = default;

  RadialGrid(RTransform *rt, double rad);

  ~RadialGrid() {
    delete[] points;
    delete[] weights;
  };

  double integrate(const int segments, double *f);

  unsigned int get_size() { return size; };
  double get_radii() { return radii; };
  double *get_points() { return points; };
  double *get_weights() { return weights; };
};

#ifdef __cplusplus
extern "C" {
#endif

RadialGrid *RadialGrid_new(RTransform *rtf, double rad);

void RadialGrid_del(RadialGrid *grid);

double RadialGrid_integrate(RadialGrid *grid, int segments, double *f);

int RadialGrid_get_size(RadialGrid *grid);

double RadialGrid_get_radii(RadialGrid *grid);

double *RadialGrid_get_points(RadialGrid *grid);

double *RadialGrid_get_weights(RadialGrid *grid);

#ifdef __cplusplus
}
#endif

#endif