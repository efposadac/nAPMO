/*file: lebedev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef ANGULAR_H
#define ANGULAR_H

#include "sphere_lebedev_rule.h"
#include "spherical_harmonics.h"
#include "utils.h"

struct _angular {
  int lorder;      // Order of quadrature.
  double *points;  // points of angular quadrature.
  double *weights; // Weights of angular quadrature
};

typedef struct _angular AngularGrid;

/*
returns the lorder number of abscissas and weights of the Lebedev quadrature
in spherical coordinates (r, theta, phi)
*/
void angular_to_spherical(AngularGrid *grid);

/*
returns the lorder number of abscissas and weights of the Lebedev quadrature in
Cartesian coordinates
*/
void angular_cartesian(AngularGrid *grid);

/*
Computes the multipolar expansion for a given function ``F``

Args:
    lmax: the order of the expansion.
    size_f: size of array f.
    f: array with the values of the function ``F`` in each point of the grid
    output: array with size :math:`lsize * size`
*/

void angular_spherical_expansion(AngularGrid *grid, const int lmax,
                                 const int size_f, double *f, double *output);

/*
Computes the multipolar expansion for a given function ``F``

Args:
    lmax: the order of the expansion.
    size_f: size of array f.
    decomposition: decomposition obtained trough angular_spherical_expansion
function.
    output: array with size :math:`lorder * size_f`

*/
void angular_eval_expansion(AngularGrid *grid, const int lmax, const int size_f,
                            double *decomposition, double *output);
/*
Integration over unit sphere :math:`4\pi` using Lebedev quadrature of ``lorder``
points.

Args:
    segments: number of appended arrays on ``f``.
    f: array with the values of the function ``F`` in each point of Lebedev's
sphere.
        w: array with the Lebedev's weights.
*/
double angular_integrate(AngularGrid *grid, const int segments, double *f);

#ifdef _CUDA

#include "cuda_helper.cuh"

/*
Convert angular quadrature from spherical to cartesian coordinates (on device)
*/
__global__ void angular_to_cartesian_cuda(const int2 gridDim, double2 *xy,
                                          double2 *zw);

#endif

#endif