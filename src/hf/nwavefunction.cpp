/*
file: nwavefunction.cpp
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#include "wavefunction.h"
#include "xc.h"

#define INDEX(i, j)                                                            \
  ((i > j) ? (((i) * ((i) + 1) / 2) + (j)) : (((j) * ((j) + 1) / 2) + (i)))

void nwavefunction_compute_2body_matrix_atm(WaveFunction *psi, BeckeGrid *grid,
                                            double *phi, double *J, double *K) {

  unsigned int ndim = psi->ndim;
  double factor = psi->x_factor;
  unsigned int size = ndim * (ndim + 1) / 2;

  // Precompute Psi
  MMap F(phi, ndim, grid->get_size());
  Matrix Psi(size, grid->get_size());

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
    }
  }

  // Compute coulomb
  A1DMap C(J, grid->get_size());
  Matrix coulomb(ndim, ndim);
  coulomb.setZero();

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= C;

      coulomb(i, j) = grid->integrate(buff);
    }
  }

  coulomb += coulomb.triangularView<Eigen::StrictlyUpper>().transpose();

  // Compute exchange
  MMap E(K, ndim, grid->get_size());
  Matrix exchange(ndim, ndim);
  exchange.setZero();

 // FELIX: add conditional for HF or hybrid functionals
  if (factor != 0.0){
    for (unsigned int i = 0; i < ndim; ++i) {
      for (unsigned int j = i; j < ndim; ++j) {

	Array1D buff = F.row(i).array() * E.row(j).array();

	exchange(i, j) = grid->integrate(buff);
      }
    }

    exchange += exchange.triangularView<Eigen::StrictlyUpper>().transpose();
    exchange *= factor;
  }

  // std::cout<<exchange<<std::endl;

  MMap G(psi->G, ndim, ndim);

  G = coulomb + exchange;
}

void nwavefunction_compute_2body_matrix_mol(WaveFunction *psi, BeckeGrid *grid,
                                            double *phi, double *J, double *K) {

  unsigned int ndim = psi->ndim;
  double factor = psi->x_factor;
  unsigned int size = ndim * (ndim + 1) / 2;

  // Precompute Psi
  MMap F(phi, ndim, grid->get_size());
  Matrix Psi(size, grid->get_size());

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
    }
  }

  // Compute coulomb
  A1DMap C(J, grid->get_size());
  Matrix coulomb(ndim, ndim);
  coulomb.setZero();

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= C;

      coulomb(i, j) = grid->integrate(buff);
    }
  }

  coulomb += coulomb.triangularView<Eigen::StrictlyUpper>().transpose();

  MMap E(K, ndim, grid->get_size());
  Matrix exchange(ndim, ndim);
  exchange.setZero();

 // FELIX: add conditional for HF or hybrid functionals
  if (factor != 0.0){
    for (unsigned int i = 0; i < ndim; ++i) {
      for (unsigned int j = i; j < ndim; ++j) {

	Array1D buff = F.row(i).array() * E.row(j).array();

	exchange(i, j) = grid->integrate(buff);
      }
    }

    exchange += exchange.triangularView<Eigen::StrictlyUpper>().transpose();
    exchange *= factor * psi->eta;
  }

  // std::cout<<exchange<<std::endl;

  MMap G(psi->G, ndim, ndim);

  G = coulomb + exchange;
}

void nwavefunction_compute_coupling(WaveFunction *psi, BeckeGrid *grid,
                                    double *phi, double *other_J, double *res) {

  unsigned int ndim = psi->ndim;
  unsigned int size = ndim * (ndim + 1) / 2;

  // Precompute <ab| in the integral <ab|r12|AB> "Coupling integral"
  MMap PHI(phi, ndim, grid->get_size());
  Matrix Psi(size, grid->get_size());

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = PHI.row(i).array() * PHI.row(j).array();
    }
  }

  // Compute coupling
  A1DMap C(other_J, grid->get_size());
  MMap coupling(res, ndim, ndim);
  coupling.setZero();

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= C;

      coupling(i, j) = grid->integrate(buff);
    }
  }

  coupling += coupling.triangularView<Eigen::StrictlyUpper>().transpose();
}

void nwavefunction_compute_density_from_dm(BasisSet *basis, BeckeGrid *grid,
                                           double *dm, double *output,
                                           double epsilon, double *dmmaxrow) {
  int nbasis = basis->get_nbasis();

  for (int point = 0; point < grid->get_size(); ++point) {

    auto aux = basis->compute(&grid->get_points()[point * 3]);

    if (epsilon > 0) {

      // compute the maximum basis function
      double absmax_basis = 0.0;
      for (int i = 0; i < nbasis; i++) {
        absmax_basis = std::max(fabs(aux[i]), absmax_basis);
      }

      // upper estimate of the density
      double rho_upper = 0.0;
      for (int i = 0; i < nbasis; i++) {
        rho_upper += fabs(aux[i]) * dmmaxrow[i];
      }
      rho_upper *= nbasis * absmax_basis;

      // if the upper bound is too low, do not compute density.
      if (rho_upper < epsilon)
        return;

      // modify epsilon to avoid recomputation
      epsilon /= absmax_basis * nbasis * nbasis;
    }

    // Loop over all basis functions and add significant contributions
    double rho = 0.0;
    for (int i = 0; i < nbasis; i++) {

      // if the contribution of this loop is smaller than
      // epsilon/nbasis, skipt
      // it.
      if (epsilon > 0) {
        if (fabs(aux[i]) * dmmaxrow[i] < epsilon)
          continue;
      }

      double tmp = 0;
      for (int j = i - 1; j >= 0; j--) {
        tmp += aux[j] * dm[i * nbasis + j];
      }
      rho += (2.0 * tmp + dm[i * (nbasis + 1)] * aux[i]) * aux[i];
    }

    output[point] = rho;
  }
}

void nwavefunction_compute_exccor_matrix(WaveFunction *psi, BeckeGrid *grid, double *phi,
                                    double *rho, double *XC) {

  unsigned int n = grid->get_size(); //grid points 
  unsigned int ndim = psi->ndim; // number of orbitals
  unsigned int size = ndim * (ndim + 1) / 2; // triangular matrix size

  // Creation of the required objects:
  // rho is density in grid
  
  // Orbitals in grid - organized
  MMap F(phi, ndim, grid->get_size());

  // Orbital products in grid - calculated
  Matrix Psi(size, grid->get_size());
  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
    }
  }

  //TODO: Density and orbitals gradients
  
  // Output: total exchange correlation energy for species
  psi->xc_energy = 0.0;
  // Output: exchange correlation matrix for species
  MMap exccor(psi->XC, ndim, ndim);
  exccor.setZero();
  // Output: exchange correlation potential grid for species
  A1DMap P(XC, grid->get_size());
  P.setZero();

  //Libxc objects
  xc_func_type exc_func;
  xc_func_type cor_func;
  int exc_func_id;
  int cor_func_id;

  // local variables to store the exchange and correlation contributions
  double *ene ; 
  ene = (double*)malloc(sizeof(double)*n); // energy density
  double *pot ; 
  pot = (double*)malloc(sizeof(double)*n); // density potential
  double *sig ;  
  sig = (double*)malloc(sizeof(double)*n); // gradient squared
  double *sigpot ; 
  sigpot = (double*)malloc(sizeof(double)*n); // gradient squared potential

  Array1D buff(n); //auxiliary for integration

  
  //Temporal
  //TODO set variables in the input 
  exc_func_id= xc_functional_get_number("LDA_X");
  cor_func_id= xc_functional_get_number("LDA_C_VWN");

  //TODO select closed shell or open shell
  xc_func_init(&exc_func, exc_func_id , XC_UNPOLARIZED);
  xc_func_init(&cor_func, cor_func_id , XC_UNPOLARIZED);

  //TODO place this somewhere else
  // printf("Grid size '%i' \n", n);
  // printf("Electron exchange functional is '%s' \n", exc_func.info->name);
  // printf("Electron exchange family is '%d' \n", exc_func.info->family);
  // printf("Electron correlation functional is '%s' \n", cor_func.info->name);
  // printf("Electron correlation family is '%d' \n", cor_func.info->family);
 
  // printf("density test '%f' \n", grid->integrate(rho));
  
  // Exchange contributions
  switch(exc_func.info->family)
    {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&exc_func, n, rho, ene, pot);
      break;
    case XC_FAMILY_GGA:
      xc_gga_exc_vxc(&exc_func, n, rho, sig, ene, pot, sigpot);
      break;
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&exc_func, n, rho, sig, ene, pot, sigpot);
      break;
    }

  for (unsigned int i = 0; i < n; ++i) {
    buff[i]=rho[i]*ene[i];
    P[i]=pot[i];
  }

  // printf("exchange energy'%f' \n", grid->integrate(buff));

  psi->xc_energy += grid->integrate(buff) ;
    
  // Correlation contributions
  memset(ene,0,n);
  memset(pot,0,n);

  switch(cor_func.info->family)
    {
    case XC_FAMILY_LDA:
      xc_lda_exc_vxc(&cor_func, n, rho, ene, pot);
      break;
    case XC_FAMILY_GGA:
      xc_gga_exc_vxc(&cor_func, n, rho, sig, ene, pot, sigpot);
      break;
    case XC_FAMILY_HYB_GGA:
      xc_gga_exc_vxc(&cor_func, n, rho, sig, ene, pot, sigpot);
      break;
    }

  for (unsigned int i = 0; i < n; ++i) {
    buff[i]=rho[i]*ene[i];
    P[i]+=pot[i];
  }

  // printf("correlation energy'%f' \n", grid->integrate(buff));

  psi->xc_energy += grid->integrate(buff) ;

  
  //Building matrix

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      buff *= P;
      
      exccor(i, j) = grid->integrate(buff);
    }
  }

  exccor += exccor.triangularView<Eigen::StrictlyUpper>().transpose();

  // tests
  // printf("xc_energy'%f' \n", psi->xc_energy) ;

  // std::cout << "\n\t exccor:\n";
  // std::cout << exccor << std::endl;
 
  // printf("XC \n");
  // for (unsigned int i = 0; i < n; ++i) {
  //   printf("%f\n", XC[i] );
  // }
 
  free(ene);
  free(pot);
  free(sig);
  free(sigpot);
  

}

// void nwavefunction_compute_exccor_matrix(WaveFunction *psi, BeckeGrid *grid, double *phi,
//                                     double *rho, double *XC) {

//   unsigned int n = grid->get_size(); //grid points 
//   unsigned int ndim = psi->ndim; // number of orbitals
//   unsigned int size = ndim * (ndim + 1) / 2; // triangular matrix size

//   // Creation of the required objects:

//   // rho is density in grid
  
//   // Orbitals in grid - organized
//   MMap F(phi, ndim, grid->get_size());

//   // Orbital products in grid - calculated
//   Matrix Psi(size, grid->get_size());
//   for (unsigned int i = 0; i < ndim; ++i) {
//     for (unsigned int j = i; j < ndim; ++j) {
//       Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
//     }
//   }

//   //TODO: Density and orbitals gradients
  
//   // Output: total exchange correlation energy for species
//   psi->xc_energy = 0.0;
//   // Output: exchange correlation matrix for species
//   MMap exccor(psi->XC, ndim, ndim);
//   exccor.setZero();
//   // Output: exchange correlation potential grid for species
//   A1DMap P(XC, grid->get_size());
//   P.setZero();

//   //Libxc objects
//   xc_func_type exc_func;
//   xc_func_type cor_func;
//   int exc_func_id;
//   int cor_func_id;

//   // local variables to store the exchange and correlation contributions
//   double *ene ; 
//   ene = (double*)malloc(sizeof(double)*n); // energy density
//   double *pot ; 
//   pot = (double*)malloc(sizeof(double)*n); // density potential
//   double *sig ;  
//   sig = (double*)malloc(sizeof(double)*n); // gradient squared
//   double *sigpot ; 
//   sigpot = (double*)malloc(sizeof(double)*n); // gradient squared potential

//   Array1D buff(n); //auxiliary for integration

  
//   //Temporal
//   //TODO set variables in the input 
//   exc_func_id= xc_functional_get_number("LDA_X");
//   cor_func_id= xc_functional_get_number("LDA_C_VWN");

//   //TODO select closed shell or open shell
//   xc_func_init(&exc_func, exc_func_id , XC_UNPOLARIZED);
//   xc_func_init(&cor_func, cor_func_id , XC_UNPOLARIZED);

//   //TODO place this somewhere else
//   // printf("Grid size '%i' \n", n);
//   // printf("Electron exchange functional is '%s' \n", exc_func.info->name);
//   // printf("Electron exchange family is '%d' \n", exc_func.info->family);
//   // printf("Electron correlation functional is '%s' \n", cor_func.info->name);
//   // printf("Electron correlation family is '%d' \n", cor_func.info->family);
 
//   // printf("density test '%f' \n", grid->integrate(rho));
  
//   // Exchange contributions
//   switch(exc_func.info->family)
//     {
//     case XC_FAMILY_LDA:
//       xc_lda_exc_vxc(&exc_func, n, rho, ene, pot);
//       break;
//     case XC_FAMILY_GGA:
//       xc_gga_exc_vxc(&exc_func, n, rho, sig, ene, pot, sigpot);
//       break;
//     case XC_FAMILY_HYB_GGA:
//       xc_gga_exc_vxc(&exc_func, n, rho, sig, ene, pot, sigpot);
//       break;
//     }

//   for (unsigned int i = 0; i < n; ++i) {
//     buff[i]=rho[i]*ene[i];
//     P[i]=pot[i];
//   }

//   // printf("exchange energy'%f' \n", grid->integrate(buff));

//   psi->xc_energy += grid->integrate(buff) ;
    
//   // Correlation contributions
//   memset(ene,0,n);
//   memset(pot,0,n);

//   switch(cor_func.info->family)
//     {
//     case XC_FAMILY_LDA:
//       xc_lda_exc_vxc(&cor_func, n, rho, ene, pot);
//       break;
//     case XC_FAMILY_GGA:
//       xc_gga_exc_vxc(&cor_func, n, rho, sig, ene, pot, sigpot);
//       break;
//     case XC_FAMILY_HYB_GGA:
//       xc_gga_exc_vxc(&cor_func, n, rho, sig, ene, pot, sigpot);
//       break;
//     }

//   for (unsigned int i = 0; i < n; ++i) {
//     buff[i]=rho[i]*ene[i];
//     P[i]+=pot[i];
//   }

//   // printf("correlation energy'%f' \n", grid->integrate(buff));

//   psi->xc_energy += grid->integrate(buff) ;

  
//   //Building matrix

//   for (unsigned int i = 0; i < ndim; ++i) {
//     for (unsigned int j = i; j < ndim; ++j) {

//       Array1D buff = Psi.row(INDEX(i, j));
//       buff *= P;
      
//       exccor(i, j) = grid->integrate(buff);
//     }
//   }

//   exccor += exccor.triangularView<Eigen::StrictlyUpper>().transpose();

//   // tests
//   // printf("xc_energy'%f' \n", psi->xc_energy) ;

//   // std::cout << "\n\t exccor:\n";
//   // std::cout << exccor << std::endl;
 
//   // printf("XC \n");
//   // for (unsigned int i = 0; i < n; ++i) {
//   //   printf("%f\n", XC[i] );
//   // }
 
//   free(ene);
//   free(pot);
//   free(sig);
//   free(sigpot);
  

// }

void nwavefunction_compute_cor2species_matrix(WaveFunction *psi, WaveFunction *otherPsi, BeckeGrid *grid, double *phi,
					      double *rho, double *otherRho, double *XC) {

  unsigned int n = grid->get_size(); //grid points 
  unsigned int ndim = psi->ndim; // number of orbitals
  unsigned int size = ndim * (ndim + 1) / 2; // triangular matrix size

  // Creation of the required objects:

  // rho is density in grid
  
  // Orbitals in grid - organized
  MMap F(phi, ndim, grid->get_size());

  // Orbital products in grid - calculated
  Matrix Psi(size, grid->get_size());
  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {
      Psi.row(INDEX(i, j)) = F.row(i).array() * F.row(j).array();
    }
  }

  //TODO: Density and orbitals gradients
  
  // Output: total exchange correlation energy for species
  // cor2speciesenergy = 0.0;
  // Output: exchange correlation matrix for species
  MMap exccor(psi->XC, ndim, ndim);
  // exccor.setZero();
  // Output: exchange correlation potential grid for species
  A1DMap P(XC, grid->get_size());

  // local variables to store the exchange and correlation contributions
  double *denominator ; 
  denominator = (double*)malloc(sizeof(double)*n); // energy denominatior
  double *ene ; 
  ene = (double*)malloc(sizeof(double)*n); // energy density
  // double *pot ; 
  // pot = (double*)malloc(sizeof(double)*n); // density potential
  Array1D pot(n); //auxiliary for integration

  
  // Energy and potential contributions

  // Parameters
  //   if(this%name .eq. "correlation:epc17-2" ) then
  // else if(this%name .eq. "correlation:epc17-1" ) then
  // STOP "The nuclear electron functional chosen is not implemented"
  double a=2.35;
  double b=2.4;
  double c=6.6;

  for (unsigned int i = 0; i < n; ++i) {
    denominator[i]=a-b*sqrt(rho[i]*otherRho[i])+c*rho[i]*otherRho[i];
    ene[i]=-rho[i]*otherRho[i]/denominator[i];
    pot[i]=(b*sqrt(rho[i])*pow(otherRho[i],3/2)-2*a*otherRho[i])/pow(denominator[i],2)/2;
  }

  
  psi->xc_energy += grid->integrate(ene)/2 ;
      
  // Potential
  for (unsigned int i = 0; i < n; ++i) {
    P[i]+=pot[i];
  }

  // Matrix
  Array2D aux(ndim,ndim);
  Array1D buff(n); //auxiliary for integration

  for (unsigned int i = 0; i < ndim; ++i) {
    for (unsigned int j = i; j < ndim; ++j) {

      Array1D buff = Psi.row(INDEX(i, j));
      
      buff *= pot;
      // for (unsigned int k = 0; i < n; ++i) {
      // 	buff[k]=Psi.row(INDEX(i, j))[k]*pot[k];
      // }
      
      aux(i,j) = grid->integrate(buff);
    }
  }

  for (unsigned int i = 0; i < ndim; ++i) {
    exccor(i,i) += aux(i,i);
    for (unsigned int j = i+1; j < ndim; ++j) {
      exccor(i,j) += aux(i,j);
      exccor(j,i) += aux(i,j);
    }
  }


}

