/*file: ints.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co
*/

#ifndef INTS_H
#define INTS_H

#include <atomic>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <sstream>
#include <stdlib.h>
#include <thread>
#include <unordered_map>
#include <vector>

/*
Buffer for integrals
*/
struct QuartetBuffer {
  std::vector<int> p;
  std::vector<int> q;
  std::vector<int> r;
  std::vector<int> s;
  std::vector<double> val;
};

#endif