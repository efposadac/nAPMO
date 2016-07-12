/*file: ints.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 0.1
efposadac@unal.edu.co
*/

#ifndef INTS_H
#define INTS_H

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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