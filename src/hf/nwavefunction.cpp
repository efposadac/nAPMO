/*
file: nwavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"

void nwavefunction_compute_density_from_dm(double* work_basis, double* dm, long nbasis, double* output, double epsilon, double* dmmaxrow) {
    
    double absmax_basis = 0.0;
    if (epsilon > 0) {
        // compute the maximum basis function
        for (long ibasis=0; ibasis<nbasis; ibasis++) {
            double tmp = fabs(work_basis[ibasis]);
            if (tmp > absmax_basis) absmax_basis = tmp;
        }
        // upper estimate of the density
        double rho_upper = 0.0;
        for (long ibasis=0; ibasis<nbasis; ibasis++) {
            rho_upper += fabs(work_basis[ibasis])*dmmaxrow[ibasis];
        }
        rho_upper *= nbasis*absmax_basis;

        // if the upper bound is too low, do not compute density.
        if (rho_upper < epsilon) return;

        // modify epsilon to avoid recomputation
        epsilon /= absmax_basis*nbasis*nbasis;
    }

    // Loop over all basis functions and add significant contributions
    double rho = 0.0;
    for (long ibasis0=0; ibasis0<nbasis; ibasis0++) {
        // if the contribution of this loop is smaller than epsilon/nbasis, skipt it.
        if (epsilon>0) {
            if (fabs(work_basis[ibasis0])*dmmaxrow[ibasis0] < epsilon)
                continue;
        }
        double tmp = 0;
        for (long ibasis1=ibasis0-1; ibasis1>=0; ibasis1--) {
            tmp += work_basis[ibasis1]*dm[ibasis0*nbasis+ibasis1];
        }
        rho += (2*tmp+dm[ibasis0*(nbasis+1)]*work_basis[ibasis0])*work_basis[ibasis0];
    }
    *output += rho;
}