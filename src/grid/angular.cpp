/*file: angular.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "angular.h"

AngularGrid::AngularGrid(const int order) : lorder(order) {
  points = new double[lorder * 3];
  weights = new double[lorder];

  cartesian();
}

void AngularGrid::cartesian() {

  double *x = new double[lorder];
  double *y = new double[lorder];
  double *z = new double[lorder];

  ld_by_order(lorder, x, y, z, weights);

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif

  for (unsigned int i = 0; i < lorder; ++i) {
    int idx = i * 3;
    points[idx + 0] = x[i];
    points[idx + 1] = y[i];
    points[idx + 2] = z[i];
  }

  delete[] x;
  delete[] y;
  delete[] z;
}

void AngularGrid::spherical() {
  double t, p;

#ifdef _OPENMP
#pragma omp parallel for default(shared) private(t, p)
#endif

  for (unsigned int i = 0; i < lorder; ++i) {
    int idx = i * 3;
    xyz_to_tp(points[idx + 0], points[idx + 1], points[idx + 2], &t, &p);

    points[idx + 0] = 1.0;
    points[idx + 1] = t;
    points[idx + 2] = p;
  }
}

void AngularGrid::spherical_expansion(const int lmax, const unsigned int size_f,
                                      double *f, double *output) {

  // Obtain Lebedev's points in spherical
  double *t = new double[lorder];
  double *p = new double[lorder];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 0; i < lorder; ++i) {
    unsigned int idx = i * 3;
    xyz_to_tp(points[idx + 0], points[idx + 1], points[idx + 2], &t[i], &p[i]);
  }

  // Calculate index for expansion.
  unsigned int lsize = (lmax + 1) * (lmax + 1);
  int *lindex = new int[lsize * 2];

  unsigned int idx = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      lindex[idx] = l;
      lindex[lsize + idx] = m;
      idx++;
    }
  }

// Calculate expansion
#ifdef _OPENMP
#pragma omp parallel default(shared)
#endif
  {
    double *s = new double[lorder];
#ifdef _OPENMP
#pragma omp for
#endif
    for (unsigned int i = 0; i < lsize; ++i) {
      spherical_harmonics_real(lindex[i], lindex[lsize + i], t, p, s, lorder);
      for (unsigned int j = 0; j < size_f; ++j) {
        output[j * lsize + i] = integrate(1, &f[j * lorder], s);
      }
    }
    delete[] s;
  }
  delete[] t;
  delete[] p;
  delete[] lindex;
}

void AngularGrid::eval_expansion(const int lmax, const unsigned int size_f,
                                 double *decomposition, double *output) {

  // Obtain Lebedev's points in spherical
  double *t = new double[lorder];
  double *p = new double[lorder];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 0; i < lorder; ++i) {
    unsigned int idx = i * 3;
    xyz_to_tp(points[idx + 0], points[idx + 1], points[idx + 2], &t[i], &p[i]);
  }

  // Calculate index for expansion.
  unsigned int lsize = (lmax + 1) * (lmax + 1);
  int *lindex = new int[lsize * 2];
  unsigned int idx = 0;
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      lindex[idx] = l;
      lindex[lsize + idx] = m;
      idx++;
    }
  }

  // Calculate spherical harmonics
  double *s = new double[lorder * lsize];

#ifdef _OPENMP
#pragma omp parallel for default(shared)
#endif
  for (unsigned int i = 0; i < lsize; ++i) {
    spherical_harmonics_real(lindex[i], lindex[lsize + i], t, p, &s[i * lorder],
                             lorder);
  }

// Evaluate expansion
#ifdef _OPENMP
#pragma omp parallel for default(shared) collapse(2)
#endif
  for (unsigned int i = 0; i < size_f; ++i) {
    for (unsigned int j = 0; j < lorder; ++j) {
      double aux = 0.0;
      for (unsigned int l = 0; l < lsize; ++l) {
        aux += decomposition[i * lsize + l] * s[l * lorder + j];
      }
      output[i * lorder + j] = aux;
    }
  }

  delete[] s;
  delete[] t;
  delete[] p;
  delete[] lindex;
}

double AngularGrid::integrate(const unsigned int segments, double *f,
                              double *s) {

  double *work;

  if (segments > 1) {
    work = new double[lorder];
    utils_multiply_segmented_array(lorder, segments, f, work);
  } else {
    work = &f[0];
  }

  double output = 0.0;
  for (unsigned int i = 0; i < lorder; ++i) {
    output += work[i] * weights[i] * s[i];
  }

  if (segments > 1) {
    delete[] work;
  }

  return output * 4.0 * M_PI;
}

/*
Python wrapper
*/

AngularGrid *AngularGrid_new(int lorder) { return new AngularGrid(lorder); }

void AngularGrid_del(AngularGrid *grid) { grid->~AngularGrid(); }

void AngularGrid_spherical(AngularGrid *grid) { grid->spherical(); }

void AngularGrid_spherical_expansion(AngularGrid *grid, const int lmax,
                                     const int size_f, double *f,
                                     double *output) {

  grid->spherical_expansion(lmax, size_f, f, output);
}

void AngularGrid_eval_expansion(AngularGrid *grid, const int lmax,
                                const int size_f, double *decomposition,
                                double *output) {

  grid->eval_expansion(lmax, size_f, decomposition, output);
}

double AngularGrid_integrate(AngularGrid *grid, const int segments, double *f,
                             double *s) {

  return grid->integrate(segments, f, s);
}

int AngularGrid_get_lorder(AngularGrid *grid) { return grid->get_lorder(); }

double *AngularGrid_get_points(AngularGrid *grid) { return grid->get_points(); }

double *AngularGrid_get_weights(AngularGrid *grid) {

  return grid->get_weights();
}