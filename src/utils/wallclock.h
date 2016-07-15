/*file: wallclock.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co*/

/*Helper function to calculate walltimes.*/

#ifndef WALLCLOCK_H
#define WALLCLOCK_H

#include <time.h>

/* compute high precision walltime and walltime difference */
inline double wallclock(const double *ref) {
  struct timespec t;
  double ret;

  clock_gettime(CLOCK_REALTIME, &t);
  ret = ((double)t.tv_sec) * 1.0e6 + 1.0e-3 * ((double)t.tv_nsec);

  return ref ? (ret - *ref) : ret;
}

#endif