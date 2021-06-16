/*
file: eigen_helper.h
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#ifndef EIGEN_HELPER
#define EIGEN_HELPER

// Eigen matrix algebra library

#define EIGEN_MATRIXBASE_PLUGIN "utils/eigen_plugins.h"

#include "utils.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/StdVector>

/*
Type definitions
*/

/*
import dense, dynamically sized Matrix type from Eigen;
this is a matrix with row-major storage
(http://en.wikipedia.org/wiki/Row-major_order)
*/
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;

/*
dynamically sized Vector type
*/
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

/*
Array 2D
*/
typedef Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic> Array2D;

/*
Array 1D
*/
typedef Eigen::Array<double, Eigen::Dynamic, 1> Array1D;
typedef Eigen::Array<long int, Eigen::Dynamic, 1> Array1Dld;

/*
Types to work with raw buffers, for vectors and matrices
*/
typedef Eigen::Map<Matrix> MMap;

typedef Eigen::Map<Vector> VMap;

typedef Eigen::Map<Array2D> A2DMap;

typedef Eigen::Map<Array1D> A1DMap;

/*
 * Compute the eigenvectors and eigenvalues, sorted
 */
template <class Derived>
template <class MATRIX1, class VECTOR1>
EIGEN_STRONG_INLINE void
Eigen::MatrixBase<Derived>::eigenVectorsVec(MATRIX1 &eVecs,
                                            VECTOR1 &eVals) const {

  // Eigen::SelfAdjointEigenSolver<Derived> es(*this);
  // eVecs = es.eigenvectors().real(); // Keep only the real part of complex
  // matrix eVals = es.eigenvalues().real();  // Keep only the real part of
  // complex matrix

  // Calculate using LAPACK
  int N = this->cols(); // Number of columns of A

  // Making A Lapack compatible
  MATRIX1 A = *this;
  A.transposeInPlace(); // for most symmetric matrices this has no effect

  double WORKDUMMY;
  int LWORK = -1; // Request optimum work size.
  int INFO = 0;

  VECTOR1 WI(N);
  MATRIX1 VL(N, N);

  // Get the optimum work size.
  char *NN = (char *)"N";
  char *VV = (char *)"V";

  dgeev_(NN, VV, &N, A.data(), &N, eVals.data(), WI.data(), VL.data(), &N,
         eVecs.data(), &N, &WORKDUMMY, &LWORK, &INFO);

  LWORK = int(WORKDUMMY);
  Vector WORK(LWORK);

  // Calculate
  dgeev_(NN, VV, &N, A.data(), &N, eVals.data(), WI.data(), VL.data(), &N,
         eVecs.data(), &N, WORK.data(), &LWORK, &INFO);

  // back to C++
  eVecs.transposeInPlace();

  // Sort by ascending eigenvalues:
  std::vector<std::pair<Scalar, Index>> D;
  D.reserve(eVals.size());

  for (Index i = 0; i < eVals.size(); i++) {
    D.push_back(std::make_pair(eVals.coeff(i, 0), i));
  }

  std::sort(D.begin(), D.end());
  MATRIX1 sortedEigs;
  sortedEigs.resizeLike(eVecs);

  for (int i = 0; i < eVals.size(); i++) {
    eVals.coeffRef(i, 0) = D[i].first;
    double m = eVecs.col(D[i].second)[0] > 0 ? 1.0 : -1.0;
    sortedEigs.col(i) = eVecs.col(D[i].second) * m;
  }

  eVecs = sortedEigs;
}

#endif