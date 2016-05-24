/*
file: libint_iface.cpp
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "include/libint_iface.h"

/*
Python functions
*/

LibintInterface *LibintInterface_new() { return new LibintInterface(); }

void LibintInterface_del(LibintInterface *lint) { lint->~LibintInterface(); }

void LibintInterface_add_particle(LibintInterface *lint, const int z,
                                  const double *center, BasisSet *basis) {
  lint->add_particle(z, center, basis);
}

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result) {

  // Default case
  libint2::Operator obtype = libint2::Operator::overlap;

  switch (integral_kind) {
  case 1: // Overlap
    obtype = libint2::Operator::overlap;
    break;
  case 2: // Kinetic
    obtype = libint2::Operator::kinetic;
    break;
  case 3: // Nuclear
    obtype = libint2::Operator::nuclear;
    break;
  }

  auto tmp = lint->compute_1body_ints(obtype);

  for (int i = 0; i < tmp.size(); ++i) {
    result[i] = tmp.array()(i);
  }
}

int LibintInterface_get_nbasis(LibintInterface *lint) {
  return lint->get_nbasis();
}

/*
LibintInterface
*/

std::vector<size_t> LibintInterface::map_shell_to_basis_function() {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell : shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

void LibintInterface::add_particle(const int z, const double *center,
                                   BasisSet *basis) {
  // add atom information
  libint2::Atom atom;
  atom.atomic_number = z;
  atom.x = center[0];
  atom.y = center[1];
  atom.z = center[2];
  atoms.push_back(atom);

  // add basis-set
  // libint2::Shell::do_enforce_unit_normalization(false);

  for (int i = 0; i < basis->n_cont; ++i) {
    int l = 0;
    for (int j = 0; j < 3; ++j) {
      l += basis->basis_l[i * 3 + j];
    }

    if (basis->basis_l[i * 3] != l)
      continue;

    int cont_size = basis->n_prim_cont[i];
    std::vector<double> exponents(cont_size);
    std::vector<double> coefficients(cont_size);

    for (int j = 0; j < cont_size; ++j) {
      int index = basis->prim_index[i];
      exponents[j] = basis->exponent[index + j];
      coefficients[j] =
          basis->coefficient[index + j];// * basis->normalization[i];
    }

    shells.push_back({{exponents},
                      {
                          {l, false, {coefficients}},
                      },
                      {{atom.x, atom.y, atom.z}}});
  }

  // nbasis, max_nprim and max_l
  nbasis = 0;
  max_nprim = 0;
  max_l = 0;
  for (const auto &shell : shells) {
    nbasis += shell.size();
    max_nprim = std::max(shell.nprim(), max_nprim);
    for (auto c : shell.contr)
      max_l = std::max(c.l, max_l);
  }
}

Matrix LibintInterface::compute_1body_ints(libint2::Operator obtype) {
  const auto n = nbasis;
  Matrix result(n, n);

  // construct the integrals engine
  libint2::Engine engine(obtype, max_nprim, max_l, 0);

  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case
  if (obtype == libint2::Operator::nuclear) {
    std::vector<std::pair<double, std::array<double, 3>>> q;
    for (const auto &atom : atoms) {
      q.push_back({static_cast<double>(atom.atomic_number),
                   {{atom.x, atom.y, atom.z}}});
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function();

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over
  // Hermitian operators: (1|2) = (2|1)
  for (unsigned int s1 = 0; s1 != shells.size(); ++s1) {
    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for (unsigned int s2 = 0; s2 <= s1; ++s2) {
      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      const auto *buf = engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding
      // blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1}
                    // block, note the transpose!
        result.block(bf2, bf1, n2, n1) = buf_mat.transpose();
    }
  }

  return result;
}
