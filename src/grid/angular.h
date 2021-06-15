/*
file: angular.h
nAPMO package
Copyright © 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#ifndef ANGULAR_H
#define ANGULAR_H

#include "../utils/spherical_harmonics.h"
#include "../utils/utils.h"

#include "sphere_lebedev_rule.h"

struct AngularGrid {

private:
  unsigned int lorder; // Order of quadrature.
  double *points;      // points of angular quadrature.
  double *weights;     // Weights of angular quadrature

  /*
  returns the lorder number of abscissas and weights of the Lebedev quadrature
  in
  Cartesian coordinates
  */
  void cartesian();

public:
  AngularGrid() : lorder(0){};

  AngularGrid(const AngularGrid &) = default;

  AngularGrid(const int order);

  ~AngularGrid() {
    delete[] points;
    delete[] weights;
  };

  /*
  returns the lorder number of abscissas and weights of the Lebedev quadrature
  in spherical coordinates (r, theta, phi)
  */
  void spherical();

  /*
  Computes the multipolar expansion for a given function ``F``

  Args:
      lmax: the order of the expansion.
      size_f: size of array f.
      f: array with the values of the function ``F`` in each point of the grid
      output: array with size :math:`lsize * size`
  */

  void spherical_expansion(const int lmax, const unsigned int size_f, double *f,
                           double *output);

  /*
  Computes the multipolar expansion for a given function ``F`` (atoms only)

  Args:
      lmax: the order of the expansion.
      size_f: size of array f.
      decomposition: decomposition obtained trough spherical_expansion
  function.
      output: array with size :math:`lorder * size_f`

  */
  void eval_expansion(const int lmax, const unsigned int size_f,
                      double *decomposition, double *output);

  /*
  Integration over unit sphere :math:`4\pi` using Lebedev quadrature of
  ``lorder``
  points.

  Args:
      segments: number of appended arrays on ``f``.
      f: array with the values of the function ``F`` in each point of Lebedev's
  sphere.
          w: array with the Lebedev's weights.
  */
  double integrate(const unsigned int segments, double *f, double *s);

  /*
  Getters
  */

  unsigned int get_lorder() { return lorder; };

  double *get_points() { return points; };

  double *get_weights() { return weights; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
AngularGrid wrapper to python
*/

AngularGrid *AngularGrid_new(int lorder);

void AngularGrid_del(AngularGrid *grid);

void AngularGrid_spherical(AngularGrid *grid);

void AngularGrid_spherical_expansion(AngularGrid *grid, const int lmax,
                                     const int size_f, double *f,
                                     double *output);

void AngularGrid_eval_expansion(AngularGrid *grid, const int lmax,
                                const int size_f, double *decomposition,
                                double *output);

double AngularGrid_integrate(AngularGrid *grid, const int segments, double *f,
                             double *s);

int AngularGrid_get_lorder(AngularGrid *grid);

double *AngularGrid_get_points(AngularGrid *grid);

double *AngularGrid_get_weights(AngularGrid *grid);

#ifdef __cplusplus
}
#endif

#ifdef _CUDA

#include "cuda_helper.cuh"

/*
Convert angular quadrature from spherical to cartesian coordinates (on device)
*/
__global__ void to_cartesian_cuda(const int2 gridDim, double2 *xy, double2 *zw);

#endif

#endif