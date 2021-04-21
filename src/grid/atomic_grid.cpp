/*file: atomic_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#include "atomic_grid.h"

AtomicGrid::AtomicGrid(AngularGrid *angular, RadialGrid *radial, double *R) {
  ang_grid = angular;
  rad_grid = radial;

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

double *AtomicGrid::spherical_expansion(double *f, int lmax) {
  /*
  Performs the spherical expansion:

  :math:`f_{\ell m}=\int_{\Omega} f(\\theta,\\varphi)\, Y_{\ell
  m}(\\theta,\\varphi)\,d\Omega`

  Args:
      lmax: Maximum :math:`\ell` order of the expansion.
      f: Array with the values of function :math:`f(x)` calculated in each point
  of the atomic grid.

  Returns:
      integral: Spherical expansion array with shape (nrad * lsize), where
  :math:`\ell_{size} = (\ell_{max} + 1)^2`
  */

  if (lmax < 0) {
    throw std::domain_error("lmax can not be negative.");
  }

  int lsize = (lmax + 1) * (lmax + 1);
  int rsize = rad_grid->get_size();
  int lorder = ang_grid->get_lorder();

  Array1Dld segments(rsize);
  segments = lorder;

  A1DMap W(weights, size);
  A1DMap F(f, size);
  Array1D FF(size);

  FF = F * W;

  double *output = new double[lsize * rsize];
  A2DMap buff(output, lsize, rsize);
  buff.setZero();

  dot_multi_moments(size, 1, FF.data(), points, origin, lmax, 4,
                    segments.data(), output, lsize);

  A1DMap sigma(integrate(1, rsize, lorder, weights), rsize);
  buff.rowwise() /= sigma.transpose();

  int counter = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      // proper norm for spherical harmonics
      buff.row(counter) *= std::sqrt(4.0 * M_PI * (2.0 * l + 1));
      counter++;
    }
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

double *AtomicGrid_spherical_expansion(AtomicGrid *grid, double *f, int lmax) {
  return grid->spherical_expansion(f, lmax);
}

int AtomicGrid_get_size(AtomicGrid *grid) { return grid->get_size(); }

double AtomicGrid_get_radii(AtomicGrid *grid) { return grid->get_radii(); }

double *AtomicGrid_get_origin(AtomicGrid *grid) { return grid->get_origin(); }

double *AtomicGrid_get_points(AtomicGrid *grid) { return grid->get_points(); }

double *AtomicGrid_get_weights(AtomicGrid *grid) { return grid->get_weights(); }
