/*file: extrapolation.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co
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

#include <cstdint>
#include <stdexcept>

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

/** @brief
        An extrapolation suitable for solutions of the Poisson equation.

    The prefactor of the left and right polynomial are such that the function
   remains
    continuous at the last grid point of the spline. The power of R on the left
   and right
    side is consistent with the boundary conditions of the solutions of the
   Poisson
    equation.
  */
class PotentialExtrapolation : public Extrapolation {
private:
  int64_t l; //!< The angular momentum for which the potential is computed.
  double amp_left;  //!< The prefactor for the polynomial for low x.
  double amp_right; //!< The prefactor for the polynomial for high x.

public:
  /** @brief
          Construct a PotentialExtrapolation.

      @param l
          The angular momentum for which the Coulomb potential is generated.
    */
  explicit PotentialExtrapolation(int64_t l);
  virtual void prepare(CubicSpline *cs); //!< Derive parameters from spline.
  virtual double eval_left(double x);    //!< Compute extrapolation for low x.
  virtual double
  eval_right(double x); //!< Derivative of extrapolation for low x.
  virtual double deriv_left(double x); //!< Compute extrapolation for high x.
  virtual double
  deriv_right(double x); //!< Derivative of extrapolation for high x.
  virtual bool has_tail() {
    return true;
  } //!< Returns true because if 1/R**(l+1) tail.

  //! The angular momentum of the Coulomb potential spline.
  int64_t get_l() { return l; }

  //! The prefactor for the polynomial for low x.
  double get_amp_left() { return amp_left; }
  //! The prefactor for the polynomial for high x.
  double get_amp_right() { return amp_right; }
};

#ifdef __cplusplus
extern "C" {
#endif

PowerExtrapolation *PowerExtrapolation_new(double power);
double PowerExtrapolation_get_power(PowerExtrapolation *extrapolation);

PotentialExtrapolation *PotentialExtrapolation_new(int l);

#ifdef __cplusplus
}
#endif

#endif