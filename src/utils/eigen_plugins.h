/*
file: eigen_plugins.h
nAPMO package
Copyright © 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#ifndef EIGEN_PLUGINS_H
#define EIGEN_PLUGINS_H

public:
// Implemented in eigen_helper.h
template <class MATRIX1, class VECTOR1>
EIGEN_STRONG_INLINE void eigenVectorsVec(MATRIX1 &eVecs, VECTOR1 &eVals) const;

#endif