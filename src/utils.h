/*file: utils.h
nAPMO package
Copyright (c) 2015, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@sissa.it*/

#ifndef UTILS_H
#define UTILS_H

#define max(a, b)                                                              \
  ({                                                                           \
    __typeof__(a) _a = (a);                                                    \
    __typeof__(b) _b = (b);                                                    \
    _a > _b ? _a : _b;                                                         \
  })

inline int utils_factorial2(const int n) {
  int output = 1;

  if (n < 1) {
    return 0;
  } else {
    for (int i = n; i > 0; i -= 2) {
      output *= i;
    }
    return output;
  }
}

#endif