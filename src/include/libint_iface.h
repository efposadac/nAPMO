/*file: libint_iface.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#ifndef LIBINT_IFACE_H
#define LIBINT_IFACE_H

#include <libint2.hpp>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <iostream>
#include <sstream>

#include "basis_set.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix; // import dense, dynamically sized Matrix type from Eigen;
            // this is a matrix with row-major storage
            // (http://en.wikipedia.org/wiki/Row-major_order)

class LibintInterface {
private:
  std::vector<libint2::Atom> atoms;
  std::vector<libint2::Shell> shells;
  std::vector<size_t> map_shell_to_basis_function();
  size_t max_nprim;
  size_t nbasis;
  int max_l;

public:
  LibintInterface() { libint2::initialize(); };
  ~LibintInterface() { libint2::finalize(); };
  void add_particle(const int z, const double *center, BasisSet *basis);
  Matrix compute_1body_ints(libint2::Operator obtype);
  size_t get_max_nprim() { return max_nprim; };
  size_t get_nbasis() { return nbasis; };
  int get_max_l() { return max_l; };
};

#ifdef __cplusplus
extern "C" {
#endif

LibintInterface *LibintInterface_new();

void LibintInterface_del(LibintInterface *lint);

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center, BasisSet *basis);
void LibintInterface_compute_1body_ints(LibintInterface *lint, int integral_kind, double *result);

int LibintInterface_get_nbasis(LibintInterface *lint);

#ifdef __cplusplus
}
#endif

#endif