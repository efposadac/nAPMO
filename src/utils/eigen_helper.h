/*file: eigen_helper.c
nAPMO package
Copyright (c) 2016, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu*/

#ifndef EIGEN_HELPER
#define EIGEN_HELPER

// Eigen matrix algebra library

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

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

/*
Types to work with raw buffers, for vectors and matrices
*/
typedef Eigen::Map<Matrix> MMap;

typedef Eigen::Map<Vector> VMap;

typedef Eigen::Map<Array2D> A2DMap;

typedef Eigen::Map<Array1D> A1DMap;

#endif