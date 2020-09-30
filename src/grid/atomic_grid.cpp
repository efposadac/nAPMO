/*file: atomic_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

#include "atomic_grid.h"

AtomicGrid::AtomicGrid(AngularGrid *angular, RadialGrid *radial, double *R) {

  size = radial->get_size() * angular->get_lorder();
  radii = radial->get_radii();

  points = new double[3 * size];
  weights = new double[size];
  origin = new double[3];

  origin[0] = R[0];
  origin[1] = R[1];
  origin[2] = R[2];

  double p[3];

  for (unsigned int i = 0; i < radial->get_size(); ++i) {
    double rp = radial->get_points()[i];
    double rw = radial->get_weights()[i];

    for (unsigned int j = 0; j < angular->get_lorder(); ++j) {
      unsigned int idx = j * 3;
      p[0] = angular->get_points()[idx + 0];
      p[1] = angular->get_points()[idx + 1];
      p[2] = angular->get_points()[idx + 2];

      double aw = angular->get_weights()[j];

      unsigned int aux = i * angular->get_lorder() + j;
      unsigned int idy = aux * 3;

      points[idy + 0] = p[0] * rp + origin[0];
      points[idy + 1] = p[1] * rp + origin[1];
      points[idy + 2] = p[2] * rp + origin[2];

      weights[aux] = aw * rw;
    }
  }
}

double *AtomicGrid::integrate(const unsigned int functions,
                              const unsigned int segments,
                              const unsigned int size, double *f) {

  int total_size = size * segments;
  double *work;

  // Convert in one array of size ``total_size``
  if (functions > 1) {
    work = new double[total_size];
    utils_multiply_segmented_array(total_size, functions, f, work);
  } else {
    work = &f[0];
  }

  // Perform reduction
  double *output = new double[segments]();
  for (unsigned int i = 0; i < segments; ++i) {
    for (unsigned int j = 0; j < size; ++j) {
      output[i] += work[i * size + j];
    }
  }

  if (functions > 1) {
    delete[] work;
  }

  return output;
}

/*
Python wrapper
*/

AtomicGrid *AtomicGrid_new(AngularGrid *angular, RadialGrid *radial,
                           double *R) {

  return new AtomicGrid(angular, radial, R);
}

void AtomicGrid_del(AtomicGrid *grid) { return grid->~AtomicGrid(); }

double *AtomicGrid_integrate(AtomicGrid *grid, int nfunc, int segments,
                             int size, double *f) {

  return grid->integrate(nfunc, segments, size, f);
}

int AtomicGrid_get_size(AtomicGrid *grid) { return grid->get_size(); }

double AtomicGrid_get_radii(AtomicGrid *grid) { return grid->get_radii(); }

double *AtomicGrid_get_origin(AtomicGrid *grid) { return grid->get_origin(); }

double *AtomicGrid_get_points(AtomicGrid *grid) { return grid->get_points(); }

double *AtomicGrid_get_weights(AtomicGrid *grid) { return grid->get_weights(); }