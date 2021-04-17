/*file: radial.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

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

double *RadialGrid::deriv2(double *f) {
  /*
  7-point finite-difference second differentiation
  (1/r)(d2/dr2)rf(r) = (2/r)(df/dr)+(d2f/dr2) =
  (2/r)(dx/dr)(df/dx) + (d2x/dr2)(df/dx) + (dx/dr)^2(d2f/dx2)
  */

  double H = M_PI / (size + 1.0);
  A1DMap FF(f, size);

  Array1D F1D(size);
  Array1D F2D(size);

  // First Derivatives
  F1D(0) =
      (-1764.0 * FF(size - 1) + 4320.0 * FF(size - 2) - 5400.0 * FF(size - 3) +
       4800.0 * FF(size - 4) - 2700.0 * FF(size - 5) + 864.0 * FF(size - 6) -
       120.0 * FF(size - 7)) /
      (720.0 * H);

  F1D(1) = (-120.0 * FF(size - 1) - 924.0 * FF(size - 2) +
            1800.0 * FF(size - 3) - 1200.0 * FF(size - 4) +
            600.0 * FF(size - 5) - 180.0 * FF(size - 6) + 24.0 * FF(size - 7)) /
           (720.0 * H);

  F1D(2) = (24.0 * FF(size - 1) - 288.0 * FF(size - 2) - 420.0 * FF(size - 3) +
            960.0 * FF(size - 4) - 360.0 * FF(size - 5) + 96.0 * FF(size - 6) -
            12.0 * FF(size - 7)) /
           (720.0 * H);

  F1D(3) = (-12.0 * FF(size - 1) + 108.0 * FF(size - 2) - 540.0 * FF(size - 3) +
            540.0 * FF(size - 5) - 108.0 * FF(size - 6) + 12.0 * FF(size - 7)) /
           (720.0 * H);

  F1D(size - 4) = (-12.0 * FF(6) + 108.0 * FF(5) - 540.0 * FF(4) +
                   540.0 * FF(2) - 108.0 * FF(1) + 12.0 * FF(0)) /
                  (720.0 * H);

  F1D(size - 3) = (12.0 * FF(6) - 96.0 * FF(5) + 360.0 * FF(4) - 960.0 * FF(3) +
                   420.0 * FF(2) + 288.0 * FF(1) - 24.0 * FF(0)) /
                  (720.0 * H);

  F1D(size - 2) =
      (-24.0 * FF(6) + 180.0 * FF(5) - 600.0 * FF(4) + 1200.0 * FF(3) -
       1800.0 * FF(2) + 924.0 * FF(1) + 120.0 * FF(0)) /
      (720.0 * H);

  F1D(size - 1) =
      (120.0 * FF(6) - 864.0 * FF(5) + 2700.0 * FF(4) - 4800.0 * FF(3) +
       5400.0 * FF(2) - 4320.0 * FF(1) + 1764.0 * FF(0)) /
      (720.0 * H);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 4; i < size - 4; ++i) {
    int idx = size - i - 1;
    F1D(i) =
        (+144.0 * FF(idx + 4) - 1536.0 * FF(idx + 3) + 8064.0 * FF(idx + 2) -
         32256.0 * FF(idx + 1) + 32256.0 * FF(idx - 1) - 8064.0 * FF(idx - 2) +
         1536.0 * FF(idx - 3) - 144.0 * FF(idx - 4)) /
        (40320.0 * H);
  }

  // Second Derivatives
  F2D(0) =
      (1624.0 * FF(size - 1) - 6264.0 * FF(size - 2) + 10530.0 * FF(size - 3) -
       10160.0 * FF(size - 4) + 5940.0 * FF(size - 5) - 1944.0 * FF(size - 6) +
       274.0 * FF(size - 7)) /
      (360.0 * H * H);

  F2D(1) = (274.0 * FF(size - 1) - 294.0 * FF(size - 2) - 510.0 * FF(size - 3) +
            940.0 * FF(size - 4) - 570.0 * FF(size - 5) + 186.0 * FF(size - 6) -
            26.0 * FF(size - 7)) /
           (360.0 * H * H);

  F2D(2) = (-26.0 * FF(size - 1) + 456.0 * FF(size - 2) - 840.0 * FF(size - 3) +
            400.0 * FF(size - 4) + 30.0 * FF(size - 5) - 24.0 * FF(size - 6) +
            4.0 * FF(size - 7)) /
           (360.0 * H * H);

  F2D(3) = (4.0 * FF(size - 1) - 54.0 * FF(size - 2) + 540.0 * FF(size - 3) -
            980.0 * FF(size - 4) + 540.0 * FF(size - 5) - 54.0 * FF(size - 6) +
            4.0 * FF(size - 7)) /
           (360.0 * H * H);

  F2D(size - 4) = (4.0 * FF(6) - 54.0 * FF(5) + 540.0 * FF(4) - 980.0 * FF(3) +
                   540.0 * FF(2) - 54.0 * FF(1) + 4.0 * FF(0)) /
                  (360.0 * H * H);

  F2D(size - 3) = (4.0 * FF(6) - 24.0 * FF(5) + 30.0 * FF(4) + 400.0 * FF(3) -
                   840.0 * FF(2) + 456.0 * FF(1) - 26.0 * FF(0)) /
                  (360.0 * H * H);

  F2D(size - 2) =
      (-26.0 * FF(6) + 186.0 * FF(5) - 570.0 * FF(4) + 940.0 * FF(3) -
       510.0 * FF(2) - 294.0 * FF(1) + 274.0 * FF(0)) /
      (360.0 * H * H);

  F2D(size - 1) =
      (274.0 * FF(6) - 1944.0 * FF(5) + 5940.0 * FF(4) - 10160.0 * FF(3) +
       10530.0 * FF(2) - 6264.0 * FF(1) + 1624.0 * FF(0)) /
      (360.0 * H * H);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 4; i < size - 4; ++i) {
    int idx = size - i - 1;
    F2D(i) =
        (-36.0 * FF(idx + 4) + 512.0 * FF(idx + 3) - 4032.0 * FF(idx + 2) +
         32256.0 * FF(idx + 1) - 57400.0 * FF(idx) + 32256.0 * FF(idx - 1) -
         4032.0 * FF(idx - 2) + 512.0 * FF(idx - 3) - 36.0 * FF(idx - 4)) /
        (20160.0 * H * H);
  }

  double radii2 = radii * radii;
  double *D = new double[size];
#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 0; i < size; ++i) {
    int idx = size - i - 1;
    double omega = (i + 1) * H;
    double X = cos(omega);
    double Y = sin(omega);
    double Y2 = Y * Y * radii2;
    double Y3 = Y2 * Y;
    double pre2 = (1.0 - X) * (1.0 - X);
    double pre3 = pre2 * (1.0 - X);
    double pre4 = pre3 * (1.0 - X);

    D[idx] =
        ((F1D(i) * ((0.5 * pre3 / (radii2 * Y)) - (0.25 * pre4 * X / Y3))) -
         (F1D(i) * pre2 / (radii * points[idx] * Y)) +
         (0.25 * F2D(i) * pre4 / Y2));
  }

  return D;
}

RadialGrid *RadialGrid_new(RTransform *rtf, double rad) {
  return new RadialGrid(rtf, rad);
}

void RadialGrid_del(RadialGrid *grid) { grid->~RadialGrid(); }

double RadialGrid_integrate(RadialGrid *grid, int segments, double *f) {
  return grid->integrate(segments, f);
}

double *RadialGrid_deriv2(RadialGrid *grid, double *f) {
  return grid->deriv2(f);
}

int RadialGrid_get_size(RadialGrid *grid) { return grid->get_size(); }

double RadialGrid_get_radii(RadialGrid *grid) { return grid->get_radii(); }

double *RadialGrid_get_points(RadialGrid *grid) { return grid->get_points(); }

double *RadialGrid_get_weights(RadialGrid *grid) { return grid->get_weights(); }