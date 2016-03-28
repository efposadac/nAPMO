/*file: extrapolation.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it
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

#ifndef EXTRAPOLATION_H
#define EXTRAPOLATION_H

class CubicSpline;

class Extrapolation {
public:
  Extrapolation(){};
  virtual ~Extrapolation(){};
  virtual void prepare(CubicSpline *cs) = 0;
  virtual double eval_left(double x) = 0;
  virtual double eval_right(double x) = 0;
  virtual double deriv_left(double x) = 0;
  virtual double deriv_right(double x) = 0;
  virtual bool has_tail() = 0;
};

#ifdef __cplusplus
extern "C" {
#endif

void Extrapolation_prepare(Extrapolation *extrapolation, CubicSpline *cs);
void Extrapolation_del(Extrapolation *extrapolation);
double Extrapolation_eval_left(Extrapolation *extrapolation, double x);
double Extrapolation_eval_right(Extrapolation *extrapolation, double x);
double Extrapolation_deriv_left(Extrapolation *extrapolation, double x);
double Extrapolation_deriv_right(Extrapolation *extrapolation, double x);
bool Extrapolation_has_tail(Extrapolation *extrapolation);

#ifdef __cplusplus
}
#endif

class ZeroExtrapolation : public Extrapolation {
public:
  virtual void prepare(CubicSpline *cs);
  virtual double eval_left(double x);
  virtual double eval_right(double x);
  virtual double deriv_left(double x);
  virtual double deriv_right(double x);
  virtual bool has_tail() { return false; };
};

#ifdef __cplusplus
extern "C" {
#endif

ZeroExtrapolation *ZeroExtrapolation_new();

#ifdef __cplusplus
}
#endif

class CuspExtrapolation : public Extrapolation {
private:
  double a0, b0, x0;

public:
  virtual void prepare(CubicSpline *cs);
  virtual double eval_left(double x);
  virtual double eval_right(double x);
  virtual double deriv_left(double x);
  virtual double deriv_right(double x);
  virtual bool has_tail() { return false; };
};

#ifdef __cplusplus
extern "C" {
#endif

CuspExtrapolation *CuspExtrapolation_new();

#ifdef __cplusplus
}
#endif

class PowerExtrapolation : public Extrapolation {
private:
  double amp, power;

public:
  PowerExtrapolation(double power) : amp(0), power(power){};
  virtual void prepare(CubicSpline *cs);
  virtual double eval_left(double x);
  virtual double eval_right(double x);
  virtual double deriv_left(double x);
  virtual double deriv_right(double x);
  virtual bool has_tail() { return true; };
  double get_power() { return power; };
};

#ifdef __cplusplus
extern "C" {
#endif

PowerExtrapolation *PowerExtrapolation_new(double power);
double PowerExtrapolation_get_power(PowerExtrapolation *extrapolation);

#ifdef __cplusplus
}
#endif

#endif