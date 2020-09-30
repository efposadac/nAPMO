/*file: ode2.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co*/

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

// UPDATELIBDOCTITLE: Second-order finite-element ODE solver using Hermite
// polynomials

#ifndef ODE2_H
#define ODE2_H

double hermite_overlap2(long xmax, long i0, bool deriv0, long i1, bool deriv1);
double hermite_overlap3(long xmax, long i0, bool deriv0, long i1, bool deriv1,
                        long i2, bool deriv2);
double hermite_node(long x, long center, bool kind, bool deriv);
double hermite_product2(long x, long i0, bool deriv0, long i1, bool deriv1);

#ifdef __cplusplus
extern "C" {
#endif

void build_ode2(double *b, double *a, double *f, double *bcs, double *coeffs,
                double *rhs, long npoint);

#ifdef __cplusplus
}
#endif

#endif
