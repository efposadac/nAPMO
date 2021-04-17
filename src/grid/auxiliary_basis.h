/*file: auxiliary_basis.h
nAPMO package
Copyright (c) 2020, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef AUXILIARY_BASIS_H
#define AUXILIARY_BASIS_H

#include <complex>
#include <iomanip>

#include "becke.h"

struct AuxiliaryBasis {
 private:
  int ncenter;
  int nao;
  int lmax;
  int halflmax;
  unsigned int max_rad;
  unsigned int max_atm;
  double lindep;
  double * aobasis;
  std::vector<int> lcenter;
  BeckeGrid * molgrid;

  void gram_schmidt(Matrix *M);

  void rylm(Array1D dest, Array1D orig, int l, int m, std::complex<double> *ylm,
            double *r);

 public:
  AuxiliaryBasis() = default;

  AuxiliaryBasis(BeckeGrid *grid, const int lm);

  ~AuxiliaryBasis() {
    ncenter = 0;
    nao = 0;
  };

  void compute_aobasis();
  double *get_aobasis() { return aobasis; };
  double *get_aobasis_center(const int center);
  unsigned int get_nao() { return nao; };
  unsigned int get_ncenter() { return ncenter; };
};

#ifdef __cplusplus
extern "C" {
#endif
AuxiliaryBasis *AuxiliaryBasis_new(BeckeGrid* grid, int lm);

void AuxiliaryBasis_del(AuxiliaryBasis *auxiliary_basis);

void AuxiliaryBasis_compute_aobasis(AuxiliaryBasis *auxiliary_basis);

double *AuxiliaryBasis_get_aobasis(AuxiliaryBasis *auxiliary_basis);

int AuxiliaryBasis_get_nao(AuxiliaryBasis *auxiliary_basis);

int AuxiliaryBasis_get_ncenter(AuxiliaryBasis *auxiliary_basis);

#ifdef __cplusplus
}
#endif

#endif