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

// UPDATELIBDOCTITLE: Evaluation of splines on grids

#ifndef EVALUATE_H
#define EVALUATE_H

#include "../utils/omp_helper.h"
#include "cell.h"
#include "cubic_spline.h"
#include "uniform.h"

#include <vector>

void eval_spline_cube(CubicSpline *spline, double *center, double *output,
                      UniformGrid *ugrid);

void eval_spline_grid(CubicSpline *spline, double *center, double *output,
                      double *points, Cell *cell, long npoint);

#ifdef __cplusplus
extern "C" {
#endif

void eval_decomposition_grid(CubicSpline **splines, double *center,
                             double *output, double *points, Cell *cell,
                             long nspline, long npoint);

#ifdef __cplusplus
}
#endif

#endif
