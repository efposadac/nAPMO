/*file: gto.cpp
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

#include "gto.h"

/*
PrimitiveGaussian Implementation
*/

PrimitiveGaussian::PrimitiveGaussian(int *ll, double *A, double z, double c)
    : zeta(z), coeff(c), norma(1.0) {

  l[0] = std::move(ll[0]);
  l[1] = std::move(ll[1]);
  l[2] = std::move(ll[2]);

  origin[0] = std::move(A[0]);
  origin[1] = std::move(A[1]);
  origin[2] = std::move(A[2]);

  normalize();
}

void PrimitiveGaussian::normalize() {
  int aux = 0;
  for (int i = 0; i < 3; ++i) {
    aux += l[i];
  }

  norma = (pow((2.0 * zeta / M_PI), 0.75)) /
          sqrt(utils_factorial2(abs(2 * l[0] - 1)) *
               utils_factorial2(abs(2 * l[1] - 1)) *
               utils_factorial2(abs(2 * l[2] - 1)) / (pow((4.0 * zeta), aux)));
}

double PrimitiveGaussian::compute(double *r) {

  double RP2 = 0.0;
  double factor = 1.0;
  for (int i = 0; i < 3; ++i) {
    double aux = r[i] - origin[i];
    factor *= pow(aux, l[i]);
    RP2 += aux * aux;
  }

  return coeff * norma * factor * exp(-zeta * RP2);
}

double PrimitiveGaussian::overlap(PrimitiveGaussian *other) {
  return overlap_primitive(this, other);
}

/*
Python wrapper
*/

PrimitiveGaussian *PrimitiveGaussian_new(int *ll, double *A, double z,
                                         double c) {
  return new PrimitiveGaussian(ll, A, z, c);
}

void PrimitiveGaussian_get_l(PrimitiveGaussian *p, int *ll) {

  ll[0] = p->get_l()[0];
  ll[1] = p->get_l()[1];
  ll[2] = p->get_l()[2];
}

void PrimitiveGaussian_get_origin(PrimitiveGaussian *p, double *A) {
  A[0] = p->get_origin()[0];
  A[1] = p->get_origin()[1];
  A[2] = p->get_origin()[2];
}

void PrimitiveGaussian_compute(PrimitiveGaussian *p, double *r, double *output,
                               int size) {
#ifdef _OMP
#pragma omp parallel for default(shared)
#endif
  for (int i = 0; i < size; ++i) {
    output[i] = p->compute(&r[i * 3]);
  }
}

double PrimitiveGaussian_overlap(PrimitiveGaussian *p, PrimitiveGaussian *op) {
  return p->overlap(op);
}

double PrimitiveGaussian_get_zeta(PrimitiveGaussian *p) {
  return p->get_zeta();
}

double PrimitiveGaussian_get_coeff(PrimitiveGaussian *p) {
  return p->get_coeff();
}

double PrimitiveGaussian_get_norma(PrimitiveGaussian *p) {
  return p->get_norma();
}

/*
ContractedGaussian class Implementation
*/

ContractedGaussian::ContractedGaussian(PrimitiveGaussian **primitives, int n)
    : nprim(n), norma(1.0) {

  origin[0] = primitives[0]->get_origin()[0];
  origin[1] = primitives[0]->get_origin()[1];
  origin[2] = primitives[0]->get_origin()[2];

  l[0] = primitives[0]->get_l()[0];
  l[1] = primitives[0]->get_l()[1];
  l[2] = primitives[0]->get_l()[2];

  for (int i = 0; i < n; ++i) {
    prim.push_back(*primitives[i]);
  }

  normalize();
}

void ContractedGaussian::normalize() { 

  norma = 1.0 / sqrt(overlap(this)); 
}

double ContractedGaussian::compute(double *r) {

  double output = 0.0;

  for (PrimitiveGaussian &p : prim) {
    output += p.compute(r);
  }

  output *= norma;

  return output;
}

double ContractedGaussian::overlap(ContractedGaussian *other) {
  double output = 0.0;

  for (PrimitiveGaussian &pa : prim) {
    for (PrimitiveGaussian &pb : other->get_prim()){
      output += pa.overlap(&pb);
    }
  }

  output *= norma * other->get_norma();

  return output;
}

/*
ContractedGaussian wrapper to Python
*/

ContractedGaussian *ContractedGaussian_new(PrimitiveGaussian **primitives,
                                           int n) {
  return new ContractedGaussian(primitives, n);
}

int ContractedGaussian_get_nprim(ContractedGaussian *c) {
  return c->get_nprim();
}

void ContractedGaussian_get_l(ContractedGaussian *c, int *ll) {
  ll[0] = c->get_l()[0];
  ll[1] = c->get_l()[1];
  ll[2] = c->get_l()[2];
}

void ContractedGaussian_get_origin(ContractedGaussian *c, double *A) {
  A[0] = c->get_origin()[0];
  A[1] = c->get_origin()[1];
  A[2] = c->get_origin()[2];
}

void ContractedGaussian_compute(ContractedGaussian *c, double *r,
                                double *output, int size) {
#ifdef _OMP
#pragma omp parallel for default(shared)
#endif
  for (int i = 0; i < size; ++i) {
    output[i] = c->compute(&r[i * 3]);
  }
}

double ContractedGaussian_overlap(ContractedGaussian *c,
                                  ContractedGaussian *oc) {
  return c->overlap(oc);
}

double ContractedGaussian_get_norma(ContractedGaussian *c) {
  return c->get_norma();
}
