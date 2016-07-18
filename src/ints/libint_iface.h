/*file:
    libint_iface.h
nAPMO package
Copyright(c) 2014, Edwin Fernando Posada
All rights reserved.
Version:
    0.1
efposadac@unal.edu.co
*/

#ifndef LIBINT_IFACE_H
#define LIBINT_IFACE_H

#include "../system/basis_set.h"
#include "../hf/wavefunction.h"
#include "ints.h"
#include "iterators.h"
#include <libint2.hpp>
#include <libint2/diis.h>

#define STACK_SIZE 30000

using shellpair_list_t = std::unordered_map<size_t, std::vector<size_t>>;

/*
Namespace extension for threads
*/
namespace libint2 {

unsigned int nthreads;

// / fires off \c nthreads instances of lambda in parallel
template <typename Lambda> void parallel_do(Lambda &lambda) {
#ifdef _OPENMP
#pragma omp parallel
  {
    auto thread_id = omp_get_thread_num();
    lambda(thread_id);
  }
#else // use C++11 threads
  std::vector<std::thread> threads;
  for (unsigned int thread_id = 0; thread_id != libint2::nthreads;
       ++thread_id) {
    if (thread_id != nthreads - 1)
      threads.push_back(std::thread(lambda, thread_id));
    else
      lambda(thread_id);
  } // threads_id
  for (unsigned int thread_id = 0; thread_id < nthreads - 1; ++thread_id)
    threads[thread_id].join();
#endif
}

} // end libint2 namespace

/*
Auxiliary functions
*/

template <libint2::Operator Kernel = libint2::Operator::coulomb>
Matrix compute_schwartz_ints(
    const std::vector<libint2::Shell> &bs1,
    const std::vector<libint2::Shell> &bs2 = std::vector<libint2::Shell>(),
    bool use_2norm = false, // use infty norm by default
    typename libint2::operator_traits<Kernel>::oper_params_type params =
        libint2::operator_traits<Kernel>::default_params());

shellpair_list_t compute_shellpair_list(
    const std::vector<libint2::Shell> &bs1,
    const std::vector<libint2::Shell> &_bs2 = std::vector<libint2::Shell>(),
    double threshold = 1e-12);

__inline__ void write_buffer(const QuartetBuffer &buffer,
                             std::ofstream &outfile);

/*
Class definition
*/

class LibintInterface {

private:
  int max_nprim;
  int nbasis;
  int max_l;
  int sID;
  shellpair_list_t obs_shellpair_list;
  std::vector<libint2::Atom> point_charges;
  std::vector<libint2::Shell> shells;
  std::vector<double> norma;
  Matrix compute_shellblock_norm(const Matrix &A);

public:
  explicit LibintInterface(const int id);

  ~LibintInterface() { libint2::finalize(); };

  void add_pointcharges(const int z, const double *center);

  void add_basis(BasisSet *basis);

  Matrix compute_1body_ints(libint2::Operator obtype);

  void init_2body_ints();

  std::vector<QuartetBuffer> *
  compute_2body_ints(const Matrix &D, const Matrix &Schwartz,
                     double precision = std::numeric_limits<double>::epsilon());

  void
  compute_2body_disk(const char *filename, const Matrix &D,
                     const Matrix &Schwartz,
                     double precision = std::numeric_limits<double>::epsilon());

  Matrix compute_2body_direct(
      const Matrix &D, const Matrix &Schwartz,
      double precision = std::numeric_limits<double>::epsilon());

  void compute_coupling_disk(
      LibintInterface &other, const char *filename,
      double precision = std::numeric_limits<double>::epsilon());

  Matrix compute_coupling_direct(
      LibintInterface &other, const Matrix &D, const bool permuted,
      double precision = std::numeric_limits<double>::epsilon());

  std::vector<size_t> map_shell_to_basis_function();

  std::vector<libint2::Shell> get_shells() { return shells; };

  double get_norma(int index) { return norma[index]; };
  size_t get_max_nprim() { return max_nprim; };
  size_t get_nbasis() { return nbasis; };
  int get_max_l() { return max_l; };
  int get_sID() { return sID; };
};

#ifdef __cplusplus
extern "C" {
#endif

LibintInterface *LibintInterface_new(int id);

void LibintInterface_del(LibintInterface *lint);

void LibintInterface_add_basis(LibintInterface *lint, BasisSet *basis);

void LibintInterface_add_pointcharges(LibintInterface *lint, const int z,
                                      const double *center);

void LibintInterface_compute_1body_ints(LibintInterface *lint,
                                        int integral_kind, double *result);

void LibintInterface_init_2body_ints(LibintInterface *lint);

std::vector<QuartetBuffer> *
LibintInterface_compute_2body_ints(LibintInterface *lint, double *dens);

void LibintInterface_compute_2body_direct(LibintInterface *lint, double *dens,
                                          double *result);

void LibintInterface_compute_2body_disk(LibintInterface *lint,
                                        const char *filename, double *dens);

void LibintInterface_compute_coupling_direct(LibintInterface *lint,
                                             LibintInterface *olint,
                                             double *dens, double *result);

void LibintInterface_compute_coupling_disk(LibintInterface *lint,
                                           LibintInterface *olint,
                                           const char *filename);
int LibintInterface_get_nbasis(LibintInterface *lint);

libint2::DIIS<Matrix> * LibintInterface_diis_new(int iter);

void LibintInterface_diis(libint2::DIIS<Matrix> * diis, WaveFunction * psi);

#ifdef __cplusplus
}
#endif

#endif
