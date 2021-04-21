/*file: becke_grid.c
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#include "becke.h"

BeckeGrid::BeckeGrid(AtomicGrid **grids, const int n) {

  ncenter = n;

  size = 0;
  for (unsigned int i = 0; i < ncenter; ++i) {
    size += grids[i]->get_size();
  }

  radii = new double[ncenter];
  origin = new double[ncenter * 3];
  points = new double[size * 3];
  weights = new double[size];

  atgrid.reserve(ncenter);

  for (unsigned int i = 0, idp = 0; i < ncenter; ++i) {
    unsigned int idx = i * 3;

    atgrid.push_back(*grids[i]);

    origin[idx + 0] = atgrid.back().get_origin()[0];
    origin[idx + 1] = atgrid.back().get_origin()[1];
    origin[idx + 2] = atgrid.back().get_origin()[2];

    radii[i] = atgrid.back().get_radii();

    for (unsigned int j = 0; j < atgrid.back().get_size(); ++j, ++idp) {
      unsigned int idy = idp * 3;
      unsigned int idj = j * 3;

      points[idy + 0] = atgrid.back().get_points()[idj + 0];
      points[idy + 1] = atgrid.back().get_points()[idj + 1];
      points[idy + 2] = atgrid.back().get_points()[idj + 2];

      weights[idp] = atgrid.back().get_weights()[j];
    }
  }

  compute_weights();
}

// Compatibility mode for Lowdin2 grids
BeckeGrid::BeckeGrid(double *p, const int sz, const int nc) {

  size = sz;
  ncenter = nc;

  points = new double[size * 3];
  weights = new double[size];
  becke_weights = new double[size];

  for (unsigned int i = 0; i < size; ++i) {
    unsigned int idx = i * 3;
    unsigned int idy = i * 4;

    points[idx + 0] = p[idy + 0];
    points[idx + 1] = p[idy + 1];
    points[idx + 2] = p[idy + 2];
    becke_weights[i] = p[idy + 3];
    weights[i] = 1.0;
  }
}

void BeckeGrid::compute_weights() {

  int asize = (ncenter * (ncenter + 1)) / 2;

  // Calculate internuclear distance
  double *R_ij = new double[asize];

  int offset = 0;
  for (unsigned int iatom = 0; iatom < ncenter; iatom++) {
    for (unsigned int jatom = 0; jatom <= iatom; jatom++) {
      R_ij[offset] = utils_distance(&origin[3 * iatom], &origin[3 * jatom]);
      offset += 1;
    }
  }

  // Calculate atomic size parameter
  double *a_ij = new double[asize];
  offset = 0;
  for (unsigned int iatom = 0; iatom < ncenter; iatom++) {
    for (unsigned int jatom = 0; jatom <= iatom; jatom++) {

      // Eq. A6
      double aux =
          (radii[iatom] - radii[jatom]) / (radii[iatom] + radii[jatom]);

      // eq. A5
      aux = aux / (aux * aux - 1.0);

      // eq. A3
      if (aux > 0.45) {
        aux = 0.45;
      } else if (aux < -0.45) {
        aux = -0.45;
      }

      a_ij[offset] = aux;
      offset += 1;
    }
  }

  // Computation of Becke weights
  becke_weights = new double[size];

  unsigned int order = 3;
  unsigned int npoint = size / ncenter;

  for (unsigned int atom = 0; atom < ncenter; ++atom) {

#ifdef _OPENMP
#pragma omp parallel for default(shared)                                       \
    firstprivate(atom, npoint) private(offset)
#endif
    for (unsigned int point = 0; point < npoint; ++point) {

      double sum = 0.0;
      double aux = 0.0;
      for (unsigned int iatom = 0; iatom < ncenter; ++iatom) {

        double p = 1.0;
        for (unsigned int jatom = 0; jatom < ncenter; ++jatom) {

          if (iatom == jatom)
            continue;

          // offset
          if (iatom < jatom) {
            offset = (jatom * (jatom + 1)) / 2 + iatom;
          } else {
            offset = (iatom * (iatom + 1)) / 2 + jatom;
          }

          // Internuclear distance (R_ij eq. 11)
          unsigned int idx = (atom * npoint + point) * 3;
          double r_i = utils_distance(&points[idx], &origin[iatom * 3]);
          double r_j = utils_distance(&points[idx], &origin[jatom * 3]);

          // \mu_ij eq. 11
          double mu_ij = (r_i - r_j) / R_ij[offset];

          // eq. A2
          double s = mu_ij + a_ij[offset] * ((1.0 - 2 * (iatom < jatom)) *
                                             (1.0 - mu_ij * mu_ij));

          // eq. 19 and 20
          for (unsigned int k = 1; k <= order; k++) {
            s = 0.5 * s * (3 - s * s);
          }

          s = 0.5 * (1.0 - s);
          p *= s;
        }

        if (iatom == atom)
          aux = p;
        sum += p;
      }
      // eq. 22
      becke_weights[atom * npoint + point] = aux / sum;
    }
  }

  delete[] R_ij;
  delete[] a_ij;
}

double BeckeGrid::integrate(double *f) {

  A1DMap BW(becke_weights, size);
  A1DMap W(weights, size);
  A1DMap F(f, size);

  return (F * BW * W).sum();
}

double BeckeGrid::integrate(Array1D &f) {

  A1DMap BW(becke_weights, size);
  A1DMap W(weights, size);

  return (f * BW * W).sum();
}

// double *BeckeGrid::poisson_solver(double *rho, int lmax) {
//   A1DMap Ri(rho, size);
//   A1DMap Wi(weights, size);

//   Array1D dens(size) = Ri * Wi;

//   int offset = 0;
//   // U = []
//   for (int center = 0; center < ncenter; ++center) {
//     Array1D p(atgrid[center].size) = dens [offset:offset + atgrid[center].size];

//     // Spherical expansion
//     if not sph_exp:
//             sph_expansion = atgrid.spherical_expansion(lmax, p)
//         else:
//             sph_expansion = sph_exp[i]

//         result = []
//         idx = 0

//         rgrid = atgrid.radial_grid
//         rtf = rgrid.rtransform
//         radii = rtf.radius_all()

//         b = napmo.CubicSpline(2 / radii, -2 / radii**2, rtf)

//         for l in range(lmax + 1):
//             for m in range(-l, l + 1):
//                 aux = np.array(sph_expansion[:, idx])

// #Build b
//                 rho = napmo.CubicSpline(aux, rtransform=rtf)

// #The approach followed here is obtained after substitution of
// #u = r * V in Eq.(21)in Becke's paper. After this transformation,
// #the boundary conditions can be implemented such that the output
// #is more accurate.
//                 fy = -4 * np.pi * rho.y
//                 fd = -4 * np.pi * rho.dx
//                 f = napmo.CubicSpline(fy, fd, rtf)

// #Derivation of boundary condition at rmax:
// #Multiply differential equation with r **l and integrate.Using
// #partial integration and the fact that V(r) = A / r * *(l + 1) for large
// #r, we find - (2l + 1) A = -4pi * int_0 ^ infty r * *2 r * *l rho(r) and so
// #V(rmax) = A / rmax * *(l + 1) = integrate(r * *l
// #rho(r)) / (2l + 1) / rmax **(l + 1)
//                 V_rmax = rgrid.integrate(
//                     rho.y * radii**l, 1
//                 ) / radii[-1]**(l + 1) / (2 * l + 1)

// #Derivation of boundary condition at rmin:
// #Same as for rmax, but multiply differential equation with r **(-l - 1)
// #and assume that V(r) = B * r * *l for small r.
//                 V_rmin = rgrid.integrate(
//                     rho.y * radii**(-l - 1), 1
//                 ) * radii[0]**(l) / (2 * l + 1)

//                 bcs = (V_rmin, None, V_rmax, None)

// #Build a
//                 a = napmo.CubicSpline(
//                     -l * (l + 1) * radii ** -2, 2 * l * (l + 1) * radii** -3, rtf)

// #Solve the differential equation
//                 v = napmo.solve_ode2(
//                     b, a, f, bcs, napmo.PotentialExtrapolation(l))

// #v = napmo.solve_ode2(
// #b, a, f, bcs, napmo.PowerExtrapolation(-l - 1))

//                 result.append(v)
//                 idx += 1

//         offset += atgrid.size
//         U.append(result)

//     return U
//   }
// }

/*
Python wrapper
*/

BeckeGrid *BeckeGrid_new(AtomicGrid **grids, int n) {

  return new BeckeGrid(grids, n);
}

BeckeGrid *BeckeGrid_from_points(double *p, int sz, int nc) {
  return new BeckeGrid(p, sz, nc);
}

void BeckeGrid_del(BeckeGrid *grid) { return grid->~BeckeGrid(); }

double BeckeGrid_integrate(BeckeGrid *grid, double *f) {
  return grid->integrate(f);
}

int BeckeGrid_get_ncenter(BeckeGrid *grid) { return grid->get_ncenter(); }

int BeckeGrid_get_size(BeckeGrid *grid) { return grid->get_size(); }

double *BeckeGrid_get_radii(BeckeGrid *grid) { return grid->get_radii(); }

double *BeckeGrid_get_points(BeckeGrid *grid) { return grid->get_points(); }

double *BeckeGrid_get_weights(BeckeGrid *grid) { return grid->get_weights(); }

double *BeckeGrid_get_origin(BeckeGrid *grid) { return grid->get_origin(); }

double *BeckeGrid_get_becke_weights(BeckeGrid *grid) {
  return grid->get_becke_weights();
}
