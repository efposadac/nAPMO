/*file: extrapolation.cpp
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

#include <cmath>

#include "cubic_spline.h"
#include "extrapolation.h"

/*
    Extrapolation Class
*/

void Extrapolation_prepare(Extrapolation *extrapolation, CubicSpline *cs) {
  extrapolation->prepare(cs);
}

void Extrapolation_del(Extrapolation *extrapolation) {
  extrapolation->~Extrapolation();
}

double Extrapolation_eval_left(Extrapolation *extrapolation, double x) {
  return extrapolation->eval_left(x);
}

double Extrapolation_eval_right(Extrapolation *extrapolation, double x) {
  return extrapolation->eval_right(x);
}

double Extrapolation_deriv_left(Extrapolation *extrapolation, double x) {
  return extrapolation->deriv_left(x);
}

double Extrapolation_deriv_right(Extrapolation *extrapolation, double x) {
  return extrapolation->deriv_right(x);
}

bool Extrapolation_has_tail(Extrapolation *extrapolation) {
  return extrapolation->has_tail();
}

/*
   ZeroExtrapolation class
*/

void ZeroExtrapolation::prepare(CubicSpline *cs) {}

double ZeroExtrapolation::eval_left(double x) { return 0.0; }

double ZeroExtrapolation::eval_right(double x) { return 0.0; }

double ZeroExtrapolation::deriv_left(double x) { return 0.0; }

double ZeroExtrapolation::deriv_right(double x) { return 0.0; }

ZeroExtrapolation *ZeroExtrapolation_new() { return new ZeroExtrapolation(); }

/*
   CuspExtrapolation class

   Only extrapolates for values smaller than lowest grid point. Extrapolation
   (of atomic densities) at large values is unreliable and waste of time.
*/

void CuspExtrapolation::prepare(CubicSpline *cs) {
  x0 = cs->get_first_x();
  if (fabs(cs->y[0]) < 0.00001) {
    // If there is no real cusp, don't care
    a0 = 0;
    b0 = 0;
  } else {
    RTransform *rtf = cs->get_rtransform();
    a0 = cs->y[0];
    b0 = cs->dt[0] / cs->y[0] / rtf->deriv(0);
  }
#ifdef DEBUG
  printf("PARS EXP EXTRAPOL a0=%f b0=%f x0=%f\n", a0, b0, x0);
#endif
}

double CuspExtrapolation::eval_left(double x) {
  return a0 * exp(b0 * (x - x0));
}

double CuspExtrapolation::eval_right(double x) { return 0.0; }

double CuspExtrapolation::deriv_left(double x) {
  return a0 * b0 * exp(b0 * (x - x0));
}

double CuspExtrapolation::deriv_right(double x) { return 0.0; }

CuspExtrapolation *CuspExtrapolation_new() { return new CuspExtrapolation(); }

/*
   PowerExtrapolation class

   This is used for potentials that obey some power law at large distances.
*/

void PowerExtrapolation::prepare(CubicSpline *cs) {
  double x = cs->get_last_x();

  amp = cs->y[cs->n - 1] * pow(x, -power);
}

double PowerExtrapolation::eval_left(double x) { return 0.0; }

double PowerExtrapolation::eval_right(double x) { return amp * pow(x, power); }

double PowerExtrapolation::deriv_left(double x) { return 0.0; }

double PowerExtrapolation::deriv_right(double x) {
  return amp * power * pow(x, power - 1);
}

PowerExtrapolation *PowerExtrapolation_new(double power) {
  return new PowerExtrapolation(power);
}

double PowerExtrapolation_get_power(PowerExtrapolation *extrapolation) {
  return extrapolation->get_power();
}

/*
   PotentialExtrapolation class

   This is used for potentials that obey specific trends at short and long
   distances,
   depending at the angular momentum (l) for which they are computed.
*/

PotentialExtrapolation::PotentialExtrapolation(int64_t l)
    : l(l), amp_left(0.0), amp_right(0.0) {
  if (l < 0) {
    throw std::domain_error("The argument l cannot be negative.");
  }
}

void PotentialExtrapolation::prepare(CubicSpline *cs) {
  amp_left = cs->y[0] / pow(cs->get_first_x(), l);
  amp_right = cs->y[cs->n - 1] * pow(cs->get_last_x(), l + 1);
}

double PotentialExtrapolation::eval_left(double x) {
  return amp_left * pow(x, l);
}

double PotentialExtrapolation::eval_right(double x) {
  return amp_right / pow(x, l + 1);
}

double PotentialExtrapolation::deriv_left(double x) {
  if (l == 0) {
    return 0.0;
  } else {
    return l * amp_left * pow(x, l - 1);
  }
}

double PotentialExtrapolation::deriv_right(double x) {
  return -(l + 1) * amp_right / pow(x, l + 2);
}

PotentialExtrapolation *PotentialExtrapolation_new(int l) {
  return new PotentialExtrapolation(l);
}
