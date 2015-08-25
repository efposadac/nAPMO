/*file: test.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef TEST_H
#define TEST_H

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../gauss_chebyshev.h"
#include "../becke_grid.h"
#include "../lebedev.h"
#include "../wallclock.h"

#ifdef _OMP
#include <omp.h>
#endif

#ifdef _CUDA
#include <cuda_runtime.h>
#endif

#endif
