/*
file: gto.cpp
nAPMO package
Copyright © 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

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

std::vector<double> PrimitiveGaussian::deriv(double *r) {
  double RP2 = 0.0;
  double factor = 1.0;
  double aux[3];
  std::vector<double> output(3);

  // Calculate RP2 and r (aux)
  for (int i = 0; i < 3; ++i) {
    aux[i] = r[i] - origin[i];
    RP2 += aux[i] * aux[i];
  }

  double exponential = exp(-zeta * RP2);

  // di = dx, dy, and dz
  for (int di = 0; di < 3; ++di) {
    l[di] += 1;
    factor = 1.0;

    for (int i = 0; i < 3; ++i) {
      factor *= pow(aux[i], l[i]);
    }

    l[di] -= 1;
    output[di] = -2.0 * zeta * factor;

    if (l[di] > 0) {
      l[di] -= 1;
      factor = 1.0;

      for (int i = 0; i < 3; ++i) {
        factor += pow(aux[i], l[i]);
      }

      l[di] += 1;
      output[di] += factor * l[di];
      output[di] *= coeff * norma * exponential;
    }
  }

  return output;
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
  // #ifdef _OMP
  // #pragma omp parallel for default(shared)
  // #endif
  for (int i = 0; i < size; ++i) {
    output[i] = p->compute(&r[i * 3]);
  }
}

void PrimitiveGaussian_deriv(PrimitiveGaussian *p, double *r, double *output,
                             int size) {
  for (int i = 0; i < size; ++i) {
    auto aux = p->deriv(&r[i * 3]);
    for (int j = 0; j < 3; ++j) {
      output[i * 3 + j] = aux[j];
    }
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

  // printf("Origin: %3.15f %3.15f %3.15f\n", origin[0], origin[1], origin[2]);

  l[0] = primitives[0]->get_l()[0];
  l[1] = primitives[0]->get_l()[1];
  l[2] = primitives[0]->get_l()[2];

  for (int i = 0; i < n; ++i) {
    prim.push_back(*primitives[i]);
  }

  normalize();
}

void ContractedGaussian::normalize() { norma = 1.0 / sqrt(overlap(this)); }

double ContractedGaussian::compute(double *r) {
  double output = 0.0;

  for (PrimitiveGaussian &p : prim) {
    output += p.compute(r);
  }

  output *= norma;

  return output;
}

std::vector<double> ContractedGaussian::deriv(double *r) {
  std::vector<double> output(3);

  for (PrimitiveGaussian &p : prim) {
    auto aux = p.deriv(r);
    for (int i = 0; i < 3; ++i) {
      output[i] += aux[i];
    }
  }

  for (int i = 0; i < 3; ++i) {
    output[i] *= norma;
  }

  return output;
}

double ContractedGaussian::overlap(ContractedGaussian *other) {
  double output = 0.0;

  for (PrimitiveGaussian &pa : prim) {
    for (PrimitiveGaussian &pb : other->get_prim()) {
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
  // #ifdef _OMP
  // #pragma omp parallel for default(shared)
  // #endif
  for (int i = 0; i < size; ++i) {
    output[i] = c->compute(&r[i * 3]);
  }
}

void ContractedGaussian_deriv(ContractedGaussian *c, double *r, double *output,
                              int size) {
  for (int i = 0; i < size; ++i) {
    auto aux = c->deriv(&r[i * 3]);
    for (int j = 0; j < 3; ++j) {
      output[i * 3 + j] = aux[j];
    }
  }
}

double ContractedGaussian_overlap(ContractedGaussian *c,
                                  ContractedGaussian *oc) {
  return c->overlap(oc);
}

double ContractedGaussian_get_norma(ContractedGaussian *c) {
  return c->get_norma();
}
