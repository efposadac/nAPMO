/*file: gto.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#ifndef GTO_H
#define GTO_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>

#include "../ints/overlap.h"
#include "../utils/utils.h"

/*
PrimitiveGaussian Class
*/
struct PrimitiveGaussian {

private:
  int l[3];
  double origin[3];
  double zeta;
  double coeff;
  double norma;

  /*
  Calculates the normalization constant of this primitive.
  */
  void normalize();

public:
  PrimitiveGaussian() = default;

  PrimitiveGaussian(const PrimitiveGaussian &) = default;

  PrimitiveGaussian(PrimitiveGaussian &&other)
      : zeta(std::move(other.zeta)), coeff(std::move(other.coeff)),
        norma(std::move(other.norma)) {
    l[0] = std::move(other.l[0]);
    l[1] = std::move(other.l[1]);
    l[2] = std::move(other.l[2]);
    origin[0] = std::move(other.origin[0]);
    origin[1] = std::move(other.origin[1]);
    origin[2] = std::move(other.origin[2]);
  };

  PrimitiveGaussian &operator=(const PrimitiveGaussian &) = default;

  PrimitiveGaussian &operator=(PrimitiveGaussian &&other) {
    zeta = std::move(other.zeta);
    coeff = std::move(other.coeff);
    norma = std::move(other.norma);
    l[0] = std::move(other.l[0]);
    l[1] = std::move(other.l[1]);
    l[2] = std::move(other.l[2]);
    origin[0] = std::move(other.origin[0]);
    origin[1] = std::move(other.origin[1]);
    origin[2] = std::move(other.origin[2]);
    return *this;
  }

  PrimitiveGaussian(int *ll, double *A, double z, double c);

  ~PrimitiveGaussian(){};

  /*
  Computes the value of the function at point ``r``.
  */
  double compute(double *r);

  /*
  Calculate the overlap integral between two primitives
  */

  double overlap(PrimitiveGaussian *other);

  /*
  Getters
  */
  int *get_l() { return l; };

  double *get_origin() { return origin; };

  double get_zeta() { return zeta; };

  double get_coeff() { return coeff; };

  double get_norma() { return norma; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
PrimitiveGaussian wrapper to python
*/

PrimitiveGaussian *PrimitiveGaussian_new(int *ll, double *A, double z,
                                         double c);

void PrimitiveGaussian_get_l(PrimitiveGaussian *p, int *ll);

void PrimitiveGaussian_get_origin(PrimitiveGaussian *p, double *A);

void PrimitiveGaussian_compute(PrimitiveGaussian *p, double *r, double *output,
                               int size);

double PrimitiveGaussian_overlap(PrimitiveGaussian *p, PrimitiveGaussian *op);

double PrimitiveGaussian_get_zeta(PrimitiveGaussian *p);

double PrimitiveGaussian_get_coeff(PrimitiveGaussian *p);

double PrimitiveGaussian_get_norma(PrimitiveGaussian *p);

#ifdef __cplusplus
}
#endif

/*
ContractedGaussian Class
*/

struct ContractedGaussian {
private:
  int l[3];
  int nprim;
  double norma;
  double origin[3];
  std::vector<PrimitiveGaussian> prim;

  /*
  Calculates the normalization constant of this contraction.
  */
  void normalize();

public:
  ContractedGaussian() = default;

  ContractedGaussian(const ContractedGaussian &) = default;

  ContractedGaussian(PrimitiveGaussian **primitives, int n);

  ~ContractedGaussian(){};

  /*
  Computes the value of the function at point ``r``.
  */
  double compute(double *r);

  /*
  Calculate the overlap integral between two contractions
  */
  double overlap(ContractedGaussian *other);

  /*
  Getters
  */
  int get_nprim() { return nprim; };

  int *get_l() { return l; };

  double *get_origin() { return origin; };

  double get_norma() { return norma; };

  std::vector<PrimitiveGaussian> &get_prim() { return prim; };
};

#ifdef __cplusplus
extern "C" {
#endif

/*
ContractedGaussian wrapper to python
*/

ContractedGaussian *ContractedGaussian_new(PrimitiveGaussian **primitives,
                                           int n);

int ContractedGaussian_get_nprim(ContractedGaussian *c);

void ContractedGaussian_get_l(ContractedGaussian *c, int *ll);

void ContractedGaussian_get_origin(ContractedGaussian *c, double *A);

void ContractedGaussian_compute(ContractedGaussian *c, double *r,
                                double *output, int size);

double ContractedGaussian_overlap(ContractedGaussian *c,
                                  ContractedGaussian *oc);

double ContractedGaussian_get_norma(ContractedGaussian *c);

#ifdef __cplusplus
}
#endif

#endif