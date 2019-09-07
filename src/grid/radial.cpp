/*file: radial.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "radial.h"
#include "stdio.h"

RadialGrid::RadialGrid(RTransform *rt, double rad) {
  radii = rad;
  rtf = rt;
  size = rt->get_npoint();
  points = new double[size];
  weights = new double[size];

  // Set the points
  double PI_4 = M_PI * 4.0;

  for (unsigned int i = 0; i < size; ++i) {
    double x = i;
    points[i] = rt->radius(x);
    weights[i] = PI_4 * rt->deriv(x) * points[i] * points[i];
  }
}

double RadialGrid::integrate(const int segments, double *f) {
  double *work;

  if (segments > 1) {
    work = new double[size];
    utils_multiply_segmented_array(size, segments, f, work);
  } else {
    work = &f[0];
  }

  double output = 0.0;
  for (unsigned int i = 0; i < size; ++i) {
    output += work[i] * weights[i];
  }

  if (segments > 1) {
    delete[] work;
  }

  return output;
}

RadialGrid *RadialGrid_new(RTransform *rtf, double rad) {
  return new RadialGrid(rtf, rad);
}

void RadialGrid_del(RadialGrid *grid) { grid->~RadialGrid(); }

double RadialGrid_integrate(RadialGrid *grid, int segments, double *f) {
  return grid->integrate(segments, f);
}

int RadialGrid_get_size(RadialGrid *grid) { return grid->get_size(); }

double RadialGrid_get_radii(RadialGrid *grid) { return grid->get_radii(); }

double *RadialGrid_get_points(RadialGrid *grid) { return grid->get_points(); }

double *RadialGrid_get_weights(RadialGrid *grid) { return grid->get_weights(); }