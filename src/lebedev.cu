/*file: lebedev.cu
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

// extern "C" {
// #include "include/lebedev.h"
// }

// #include "cuda_helper.cuh"

// #define THREADS_PER_BLOCK 64

// __global__ void lebedev_to_cartesian_cuda(const int2 gridDim, double2 *xy,
//                                           double2 *zw) {
//   __shared__ double2 aux_xy[THREADS_PER_BLOCK];
//   __shared__ double2 aux_zw[THREADS_PER_BLOCK];

//   double sin_t, sin_p, cos_t, cos_p;

//   const unsigned int j = __umul24(blockIdx.x, blockDim.x) + threadIdx.x;

//   if (j < gridDim.y) {
//     aux_xy[threadIdx.x] = xy[j];

//     sincos(aux_xy[threadIdx.x].x, &sin_t, &cos_t);
//     sincos(aux_xy[threadIdx.x].y, &sin_p, &cos_p);

//     aux_xy[threadIdx.x] = make_double2(sin_t * cos_p, sin_t * sin_p);
//     xy[j] = aux_xy[threadIdx.x];

//     aux_zw[threadIdx.x] = zw[j];
//     aux_zw[threadIdx.x] = make_double2(cos_t, aux_zw[threadIdx.x].y);
//     zw[j] = aux_zw[threadIdx.x];
//   }
// }
