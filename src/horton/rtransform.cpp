/*file: rtransform.cpp
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
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

//#define DEBUG

#ifdef DEBUG
#include <cstdio>
#endif

#include <cmath>
#include <stdexcept>

#include "cubic_spline.h"
#include "rtransform.h"

/*
   RTransform
*/

RTransform::RTransform(int npoint) : npoint(npoint) {
  if (npoint < 2)
    throw std::domain_error("A radial grid consists of at least two points.");
}

void RTransform::radius_array(double *t, double *r, int n) {
  while (n > 0) {
    *r = radius(*t);
    n--;
    t++;
    r++;
  }
}

void RTransform::deriv_array(double *t, double *d, int n) {
  while (n > 0) {
    *d = deriv(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::deriv2_array(double *t, double *d, int n) {
  while (n > 0) {
    *d = deriv2(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::deriv3_array(double *t, double *d, int n) {
  while (n > 0) {
    *d = deriv3(*t);
    n--;
    t++;
    d++;
  }
}

void RTransform::inv_array(double *r, double *t, int n) {
  while (n > 0) {
    *t = inv(*r);
    n--;
    r++;
    t++;
  }
}

void RTransform_del(RTransform *rtransform) { rtransform->~RTransform(); }

double RTransform_radius(RTransform *rtransform, double t) {
  return rtransform->radius(t);
}

double RTransform_deriv(RTransform *rtransform, double t) {
  return rtransform->deriv(t);
}

double RTransform_deriv2(RTransform *rtransform, double t) {
  return rtransform->deriv2(t);
}

double RTransform_deriv3(RTransform *rtransform, double t) {
  return rtransform->deriv3(t);
}

double RTransform_inv(RTransform *rtransform, double r) {
  return rtransform->inv(r);
}

void RTransform_radius_array(RTransform *rtransform, double *t, double *r,
                             int n) {
  rtransform->radius_array(t, r, n);
}

void RTransform_deriv_array(RTransform *rtransform, double *t, double *d,
                            int n) {
  rtransform->deriv_array(t, d, n);
}

void RTransform_deriv2_array(RTransform *rtransform, double *t, double *d,
                             int n) {
  rtransform->deriv2_array(t, d, n);
}

void RTransform_deriv3_array(RTransform *rtransform, double *t, double *d,
                             int n) {
  rtransform->deriv3_array(t, d, n);
}

void RTransform_inv_array(RTransform *rtransform, double *r, double *t, int n) {
  rtransform->inv_array(r, t, n);
}

int RTransform_get_npoint(RTransform *rtransform) {
  return rtransform->get_npoint();
}

/*
   IdentityRTransform
*/

double IdentityRTransform::radius(double t) { return t; }

double IdentityRTransform::deriv(double t) { return 1.0; }

double IdentityRTransform::deriv2(double t) { return 0.0; }

double IdentityRTransform::deriv3(double t) { return 0.0; }

double IdentityRTransform::inv(double r) { return r; }

IdentityRTransform *IdentityRTransform_new(int npoint) {
  return new IdentityRTransform(npoint);
}

/*
   LinearRTransform
*/

LinearRTransform::LinearRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  alpha = (rmax - rmin) / (npoint - 1);
}

double LinearRTransform::radius(double t) { return alpha * t + rmin; }

double LinearRTransform::deriv(double t) { return alpha; }

double LinearRTransform::deriv2(double t) { return 0.0; }

double LinearRTransform::deriv3(double t) { return 0.0; }

double LinearRTransform::inv(double r) { return (r - rmin) / alpha; }

LinearRTransform *LinearRTransform_new(double rmin, double rmax, int npoint) {
  return new LinearRTransform(rmin, rmax, npoint);
}

/*
   ExpRTransform
*/

ExpRTransform::ExpRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  if ((rmin <= 0.0) || (rmax <= 0.0))
    throw std::domain_error("The minimum and maximum radii must be positive.");
  alpha = log(rmax / rmin) / (npoint - 1);
}

double ExpRTransform::radius(double t) { return rmin * exp(t * alpha); }

double ExpRTransform::deriv(double t) { return rmin * alpha * exp(t * alpha); }

double ExpRTransform::deriv2(double t) {
  return rmin * alpha * alpha * exp(t * alpha);
}

double ExpRTransform::deriv3(double t) {
  return rmin * alpha * alpha * alpha * exp(t * alpha);
}

double ExpRTransform::inv(double r) { return log(r / rmin) / alpha; }

ExpRTransform *ExpRTransform_new(double rmin, double rmax, int npoint) {
  return new ExpRTransform(rmin, rmax, npoint);
}

/*
   ShiftedExpRTransform
*/

ShiftedExpRTransform::ShiftedExpRTransform(double rmin, double rshift,
                                           double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rshift(rshift), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  if ((rmin <= 0.0) || (rmax <= 0.0))
    throw std::domain_error("The minimum and maximum radii must be positive.");
  r0 = rmin + rshift;
  if (r0 <= 0.0)
    throw std::domain_error("The parameter r0 must be positive.");
  alpha = log((rmax + rshift) / r0) / (npoint - 1);
}

double ShiftedExpRTransform::radius(double t) {
  return r0 * exp(t * alpha) - rshift;
}

double ShiftedExpRTransform::deriv(double t) {
  return r0 * alpha * exp(t * alpha);
}

double ShiftedExpRTransform::deriv2(double t) {
  return r0 * alpha * alpha * exp(t * alpha);
}

double ShiftedExpRTransform::deriv3(double t) {
  return r0 * alpha * alpha * alpha * exp(t * alpha);
}

double ShiftedExpRTransform::inv(double r) {
  return log((r + rshift) / r0) / alpha;
}

ShiftedExpRTransform *ShiftedExpRTransform_new(double rmin, double rshift,
                                               double rmax, int npoint) {
  return new ShiftedExpRTransform(rmin, rshift, rmax, npoint);
}

/*
   PowerRTransform
*/

PowerRTransform::PowerRTransform(double rmin, double rmax, int npoint)
    : RTransform(npoint), rmin(rmin), rmax(rmax) {
  if (rmin >= rmax)
    throw std::domain_error("rmin must be below rmax.");
  if ((rmin <= 0.0) || (rmax <= 0.0))
    throw std::domain_error("The minimum and maximum radii must be positive.");
  power = (log(rmax) - log(rmin)) / log(npoint);
  if (power < 2.0)
    throw std::domain_error(
        "Power must be at least two for a decent intgration.");
}

double PowerRTransform::radius(double t) { return rmin * pow(t + 1, power); }

double PowerRTransform::deriv(double t) {
  return power * rmin * pow(t + 1, power - 1);
}

double PowerRTransform::deriv2(double t) {
  return power * (power - 1) * rmin * pow(t + 1, power - 2);
}

double PowerRTransform::deriv3(double t) {
  return power * (power - 1) * (power - 2) * rmin * pow(t + 1, power - 3);
}

double PowerRTransform::inv(double r) { return pow(r / rmin, 1.0 / power) - 1; }

PowerRTransform *PowerRTransform_new(double rmin, double rmax, int npoint) {
  return new PowerRTransform(rmin, rmax, npoint);
}

double PowerRTransform_get_rmin(PowerRTransform *rtransform) {
  return rtransform->get_rmin();
}

double PowerRTransform_get_rmax(PowerRTransform *rtransform) {
  return rtransform->get_rmax();
}

double PowerRTransform_get_power(PowerRTransform *rtransform) {
  return rtransform->get_power();
}

/*
   ChebyshevRTransform

   r = r_m (1 + x) / (1 - x)
   x = cos(\pi * i / (N + 1))
*/

ChebyshevRTransform::ChebyshevRTransform(double radii, int npoint)
    : RTransform(npoint), radii(radii) {
  z = M_PI / (npoint + 1);
}

double ChebyshevRTransform::radius(double t) {
  t = npoint - t;
  double aux = cos(t * z);
  return radii * (1 + aux) / (1 - aux);
}

double ChebyshevRTransform::deriv(double t) {
  t = npoint - t;
  double aux = t * z;
  return (2 * radii * z * sin(aux)) / pow(cos(aux) - 1, 2);
}

double ChebyshevRTransform::deriv2(double t) {
  t = npoint - t;
  double aux = t * z;
  return 0.5 * z * z * radii * (2 + cos(aux)) * pow(1 / sin(aux * 0.5), 4);
}

double ChebyshevRTransform::deriv3(double t) {
  t = npoint - t;
  double aux = t * z;
  double aux_2 = aux * 0.5;
  return 0.25 * z * z * z * radii * (11 * cos(aux_2) + cos(1.5 * aux)) * pow(1 / sin(aux_2), 5);
}

double ChebyshevRTransform::inv(double r) {
  return -(acos((r - radii) / (r + radii)) / z)  + npoint;
}

ChebyshevRTransform *ChebyshevRTransform_new(double radii, int npoint) {
  return new ChebyshevRTransform(radii, npoint);
}

double ChebyshevRTransform_get_radii(ChebyshevRTransform *rtransform) {
  return rtransform->get_radii();
}

double ChebyshevRTransform_get_z(ChebyshevRTransform *rtransform) {
  return rtransform->get_z();
}