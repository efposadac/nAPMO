/*file: rtransform.h
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

// Transformation from uniform 1D to non-uniform 1D grids

#ifndef RTRANSFORM_H
#define RTRANSFORM_H

class RTransform {
protected:
  int npoint;

public:
  RTransform(int npoint);
  virtual ~RTransform(){};
  virtual double radius(double t) = 0;
  virtual double deriv(double t) = 0;
  virtual double deriv2(double t) = 0;
  virtual double deriv3(double t) = 0;
  virtual double inv(double r) = 0;

  void radius_array(double *t, double *r, int n);
  void deriv_array(double *t, double *d, int n);
  void deriv2_array(double *t, double *d, int n);
  void deriv3_array(double *t, double *d, int n);
  void inv_array(double *r, double *t, int n);
  int get_npoint() { return npoint; };
};

#ifdef __cplusplus
extern "C" {
#endif

void RTransform_del(RTransform *rtransform);
double RTransform_radius(RTransform *rtransform, double t);
double RTransform_deriv(RTransform *rtransform, double t);
double RTransform_deriv2(RTransform *rtransform, double t);
double RTransform_deriv3(RTransform *rtransform, double t);
double RTransform_inv(RTransform *rtransform, double r);

void RTransform_radius_array(RTransform *rtransform, double *t, double *r,
                             int n);
void RTransform_deriv_array(RTransform *rtransform, double *t, double *d,
                            int n);
void RTransform_deriv2_array(RTransform *rtransform, double *t, double *d,
                             int n);
void RTransform_deriv3_array(RTransform *rtransform, double *t, double *d,
                             int n);
void RTransform_inv_array(RTransform *rtransform, double *r, double *t, int n);
int RTransform_get_npoint(RTransform *rtransform);

#ifdef __cplusplus
}
#endif

class IdentityRTransform : public RTransform {
public:
  IdentityRTransform(int npoint) : RTransform(npoint){};
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);
};

#ifdef __cplusplus
extern "C" {
#endif

IdentityRTransform *IdentityRTransform_new(int npoint);

#ifdef __cplusplus
}
#endif

class LinearRTransform : public RTransform {
private:
  double rmin, rmax, alpha;

public:
  LinearRTransform(double rmin, double rmax, int npoint);
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);

  double get_rmin() { return rmin; };
  double get_rmax() { return rmax; };
  double get_alpha() { return alpha; };
};

#ifdef __cplusplus
extern "C" {
#endif

LinearRTransform *LinearRTransform_new(double rmin, double rmax, int npoint);
double LinearRTransform_get_rmin(LinearRTransform *rtransform);
double LinearRTransform_get_rmax(LinearRTransform *rtransform);
double LinearRTransform_get_alpha(LinearRTransform *rtransform);

#ifdef __cplusplus
}
#endif

class ExpRTransform : public RTransform {
private:
  double rmin, rmax, alpha;

public:
  ExpRTransform(double rmin, double rmax, int npoint);
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);

  double get_rmin() { return rmin; };
  double get_rmax() { return rmax; };
  double get_alpha() { return alpha; };
};

#ifdef __cplusplus
extern "C" {
#endif

ExpRTransform *ExpRTransform_new(double rmin, double rmax, int npoint);
double ExpRTransform_get_rmin(ExpRTransform *rtransform);
double ExpRTransform_get_rmax(ExpRTransform *rtransform);
double ExpRTransform_get_alpha(ExpRTransform *rtransform);

#ifdef __cplusplus
}
#endif

class ShiftedExpRTransform : public RTransform {
private:
  double rmin, rshift, rmax, r0, alpha;

public:
  ShiftedExpRTransform(double rmin, double rshift, double rmax, int npoint);
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);

  double get_rmin() { return rmin; };
  double get_rshift() { return rshift; };
  double get_rmax() { return rmax; };
  double get_r0() { return r0; };
  double get_alpha() { return alpha; };
};

#ifdef __cplusplus
extern "C" {
#endif

ShiftedExpRTransform *ShiftedExpRTransform_new(double rmin, double rshift,
                                               double rmax, int npoint);
double ShiftedExpRTransform_get_rmin(ShiftedExpRTransform *rtransform);
double ShiftedExpRTransform_get_rshift(ShiftedExpRTransform *rtransform);
double ShiftedExpRTransform_get_rmax(ShiftedExpRTransform *rtransform);
double ShiftedExpRTransform_get_r0(ShiftedExpRTransform *rtransform);
double ShiftedExpRTransform_get_alpha(ShiftedExpRTransform *rtransform);

#ifdef __cplusplus
}
#endif

class PowerRTransform : public RTransform {
private:
  double rmin, rmax;
  double power;

public:
  PowerRTransform(double rmin, double rmax, int npoint);
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);

  double get_rmin() { return rmin; };
  double get_rmax() { return rmax; };
  double get_power() { return power; };
};

#ifdef __cplusplus
extern "C" {
#endif

PowerRTransform *PowerRTransform_new(double rmin, double rmax, int npoint);
double PowerRTransform_get_rmin(PowerRTransform *rtransform);
double PowerRTransform_get_rmax(PowerRTransform *rtransform);
double PowerRTransform_get_power(PowerRTransform *rtransform);

#ifdef __cplusplus
}
#endif

class ChebyshevRTransform : public RTransform {
private:
  double radii;
  double z;

public:
  ChebyshevRTransform(double radii, int npoint);
  virtual double radius(double t);
  virtual double deriv(double t);
  virtual double deriv2(double t);
  virtual double deriv3(double t);
  virtual double inv(double r);

  double get_radii() { return radii; };
  double get_z() { return z; };
};

#ifdef __cplusplus
extern "C" {
#endif

ChebyshevRTransform *ChebyshevRTransform_new(double radii, int npoint);
double ChebyshevRTransform_get_radii(ChebyshevRTransform *rtransform);
double ChebyshevRTransform_get_z(ChebyshevRTransform *rtransform);

#ifdef __cplusplus
}
#endif

#endif
