/*file: dot_multi.h
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2016 The HORTON Development Team
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

// UPDATELIBDOCTITLE: Auxiliary functions

#ifndef DOT_MULTI_H
#define DOT_MULTI_H

#include "uniform.h"

#ifdef __cplusplus
extern "C" {
#endif

void dot_multi(long npoint, long nvector, double **data, long *segments,
               double *output);
void dot_multi_moments_cube(long nvector, double **data, UniformGrid *ugrid,
                            double *center, long lmax, long mtype,
                            double *output, long nmoment);
void dot_multi_moments(long npoint, long nvector, double **data, double *points,
                       double *center, long lmax, long mtype, long *segments,
                       double *output, long nmoment);

#ifdef __cplusplus
}
#endif

#endif
