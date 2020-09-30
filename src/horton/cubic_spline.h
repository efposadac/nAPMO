/*file: cubic_spline.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu
*/

// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2015 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

// One-dimensional cubic splines (on uniform grids)

#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "extrapolation.h"
#include "rtransform.h"

#ifdef __cplusplus
extern "C" {
#endif

void tridiagsym_solve(double *diag_mid, double *diag_up, double *right,
                      double *solution, int n);
void solve_cubic_spline_system(double *y, double *d, int npoint);
void compute_cubic_spline_int_weights(double *weights, int npoint);

#ifdef __cplusplus
}
#endif

class Extrapolation;

class CubicSpline {
private:
  Extrapolation *extrapolation;
  RTransform *rtf;
  double first_x, last_x;

public:
  double *y;
  double *dt;
  int n;

  CubicSpline(double *y, double *dt, Extrapolation *extrapolation,
              RTransform *rtf, int n);
  void eval(const double* new_x, double* new_y, int new_n);
  void eval_deriv(double *new_x, double *new_dx, int new_n);
  RTransform *get_rtransform() { return rtf; }
  // position of first (transformed) grid point
  double get_first_x() { return first_x; };
  // position of first (transformed) last point
  double get_last_x() { return last_x; };
  Extrapolation *get_extrapolation() { return extrapolation; };
};

#ifdef __cplusplus
extern "C" {
#endif

CubicSpline *CubicSpline_new(double *y, double *dt,
                             Extrapolation *extrapolation, RTransform *rtf,
                             int n);
void CubicSpline_del(CubicSpline *cspline);
void CubicSpline_eval(CubicSpline *cspline, double *new_x, double *new_y,
                      int new_n);
void CubicSpline_eval_deriv(CubicSpline *cspline, double *new_x, double *new_dx,
                            int new_n);
double CubicSpline_get_first_x(CubicSpline *cspline);
double CubicSpline_get_last_x(CubicSpline *cspline);
RTransform *CubicSpline_get_rtransform(CubicSpline *cspline);
Extrapolation *CubicSpline_get_extrapolation(CubicSpline *cspline);

#ifdef __cplusplus
}
#endif

#endif
