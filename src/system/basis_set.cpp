/*
file: basis_set.cpp
nAPMO package
Copyright © 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#include "basis_set.h"

BasisSet::BasisSet(ContractedGaussian **contractions, int n)
    : nbasis(n), max_l(0), max_nprim(0) {

  for (int i = 0; i < n; ++i) {
    cont.push_back(*contractions[i]);
  }

  for (ContractedGaussian &c : cont) {
    max_nprim = std::max(c.get_nprim(), max_nprim);
    max_l = std::max(c.get_l()[0], max_l);
  }
}

std::vector<double> BasisSet::compute(double *r) {

  std::vector<double> res;

  for (ContractedGaussian &c : cont) {
    res.push_back(c.compute(r));
  }

  return res;
}

std::vector<double> BasisSet::deriv(double *r) {

  std::vector<double> res;

  for (ContractedGaussian &c : cont) {
    res.push_back(c.compute(r));
  }

  return res;
}

void BasisSet::update(BasisSet *other) {

  nbasis += other->get_nbasis();
  max_nprim = std::max(max_nprim, other->get_max_nprim());
  max_l = std::max(max_l, other->get_max_l());

  for (int i = 0; i < other->get_nbasis(); ++i) {
    cont.push_back(other->get_cont()[i]);
  }
}

/*
BasisSet wrapper to python
*/

BasisSet *BasisSet_new_empty() { return new BasisSet(); }

BasisSet *BasisSet_new(ContractedGaussian **primitives, int n) {

  return new BasisSet(primitives, n);
}

void BasisSet_compute(BasisSet *basis, double *r, double *output, int size) {

  int nbasis = basis->get_nbasis();

#ifdef _OMP
#pragma omp parallel for default(shared)
#endif
  for (int i = 0; i < size; ++i) {
    auto aux = basis->compute(&r[i * 3]);
    for (int j = 0; j < nbasis; ++j) {
      output[i * nbasis + j] = aux[j];
    }
  }
}

void BasisSet_deriv(BasisSet *basis, double *r, double *output, int size) {

  int nbasis = basis->get_nbasis();

#ifdef _OMP
#pragma omp parallel for default(shared)
#endif
  for (int i = 0; i < size; ++i) {
    auto aux = basis->compute(&r[i * 3]);
    for (int j = 0; j < nbasis; ++j) {
      for (int k = 0; k < 3; ++k) {
        output[(i * nbasis * 3) + (j * 3) + k] = aux[j * 3 + k];
      }
    }
  }
}

void BasisSet_update(BasisSet *basis, BasisSet *other) { basis->update(other); }

int BasisSet_get_nbasis(BasisSet *basis) { return basis->get_nbasis(); }

int BasisSet_get_max_l(BasisSet *basis) { return basis->get_max_l(); }

int BasisSet_get_max_nprim(BasisSet *basis) { return basis->get_max_nprim(); }
