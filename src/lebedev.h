/*file: lebedev.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef LEBEDEV_H
#define LEBEDEV_H

#include "sphere_lebedev_rule.h"
#include "spherical_harmonics.h"
#include "utils.h"

/*
returns the lorder number of abscissas and weights of the Lebedev quadrature in
spherical coordinates
*/
void lebedev_spherical(const int lorder, double *t, double *p, double *w);

/*
returns the lorder number of abscissas and weights of the Lebedev quadrature in
Cartesian coordinates
*/
void lebedev_cartesian(const int lorder, double *points, double *w);

/*
Computes the multipolar expansion for a given function ``F``

Args:
    lorder: order of Lebedev's quadrature.
    lmax: the order of the expansion.
    size: number of radial points.
    f: array with the values of the function ``F`` in each point of the grid
    output: array with size :math:`lsize * size`
*/

void lebedev_spherical_expansion(const int lorder, const int lmax,
                                 const int size, double *f, double *output);

/*
Integration over unit sphere :math:`4\pi` using Lebedev quadrature of ``lorder``
points.

Args:
    lorder: order of Lebedev's quadrature.
    nfunc: number of appended arrays on ``f``.
    f: array with the values of the function ``F`` in each point of Lebedev's
sphere.
        w: array with the Lebedev's weights.
*/
double lebedev_integrate(const int lorder, const int nfunc, double *f,
                         double *w);

#ifdef _CUDA

#include "cuda_helper.cuh"

/*
Convert angular quadrature from spherical to cartesian coordinates (on device)
*/
__global__ void lebedev_to_cartesian_cuda(const int2 gridDim, double2 *xy,
                                          double2 *zw);

#endif

#endif