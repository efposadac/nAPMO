/*
file: wavefunction.h
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include "../grid/becke.h"
#include "../ints/ints.h"
#include "../system/basis_set.h"
#include "../utils/eigen_helper.h"
#include "../utils/omp_helper.h"
#include <math.h>

struct wf {
  double *S;  // Overlap matrix
  double *X;  // Tranformation matrix
  double *T;  // Kinetic matrix
  double *V;  // Nuclear matrix
  double *H;  // One-particle Hamiltonian
  double *C;  // Coefficients Matrix
  double *D;  // Density matrix
  double *L;  // Last Density matrix
  double *G;  // two-particle matrix
  double *J;  // Coupling matrix
  double *XC; // Exchange-Correlation matrix
  double *F;  // Fock matrix
  double *O;  // Orbitals energy
  // Convenience Variables
  int nbasis; // number of basis
  int ndim;   // Size of the matrix
  int occupation;
  double eta;       // Constant of coupling
  double kappa;     // Constant of coupling
  double x_factor;  // Fraction of exchange
  double xc_energy; // Exchange-Correlation Energy
  double energy;    // HF Energy
  double rmsd;      // Density root-mean-square deviation
};
typedef struct wf WaveFunction;

#ifdef __cplusplus
extern "C" {
#endif

void wavefunction_guess_hcore(WaveFunction *psi);

void wavefunction_transformation_matrix(WaveFunction *psi);

void wavefunction_compute_coefficients(WaveFunction *psi);

void wavefunction_compute_density(WaveFunction *psi);

void wavefunction_compute_energy(WaveFunction *psi);

void wavefunction_iterate(WaveFunction *psi);

void wavefunction_compute_energy(WaveFunction *psi);

void wavefunction_compute_2body_matrix(WaveFunction *psi,
                                       std::vector<QuartetBuffer> *ints);

// nwavefunction

void nwavefunction_compute_density_from_dm(BasisSet *basis, BeckeGrid *grid,
                                           double *dm, double *output,
                                           double epsilon, double *dmmaxrow);

void nwavefunction_compute_2body_matrix(WaveFunction *psi, BeckeGrid *grid,
                                        double *phi, double *J, double *K);

void nwavefunction_compute_coupling(WaveFunction *psi, BeckeGrid *grid,
                                    double *phi, double *other_J, double *res);

void nwavefunction_compute_xc_matrix(WaveFunction *psi, BeckeGrid *grid,
                                     double *phi, double *xc_vrho);

#ifdef __cplusplus
}
#endif

#endif
