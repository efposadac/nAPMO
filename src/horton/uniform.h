/*file: uniform.h
nAPMO package
Copyright (c) 2014, Edwin Fernando Posada
All rights reserved.
Version: 1.0
fernando.posada@temple.edu.co
*/

// HORTON: Helpful Open-source Research TOol for N-fermion systems.
// Copyright (C) 2011-2015 The HORTON Development Team
//
// This file is part of HORTON.
//
// HORTON is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// HORTON is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>
//
//--

// UPDATELIBDOCTITLE: Uniform 3D grids

#ifndef UNIFORM_H
#define UNIFORM_H

#include "cell.h"
#include "cubic_spline.h"

class UniformGrid {
public:
  double origin[3];
  double grid_rvecs[9];
  long shape[3];
  long pbc[3];

  UniformGrid(double *_origin, double *_grid_rvecs, long *_shape, long *_pbc);

  Cell *get_cell();
  Cell *get_grid_cell();

  void set_ranges_rcut(double *center, double rcut, long *ranges_begin,
                       long *ranges_end);
  double dist_grid_point(double *center, long *i);
  void delta_grid_point(double *center, long *i);
  double *get_pointer(double *array, long *i);
};

class UniformGridWindow {
private:
  UniformGrid *ugrid;

public:
  long begin[3];
  long end[3];

  UniformGridWindow(UniformGrid *ugrid, long *_begin, long *_end);

  double *get_pointer(double *array, long *j);

  void extend(double *cell, double *local);
  void wrap(double *local, double *cell);
};

long index_wrap(long i, long high);

class Range3Iterator {
private:
  const long *ranges_begin;
  const long *ranges_end;
  const long *shape;
  long loop_shape[3];
  long npoint;

public:
  Range3Iterator(const long *ranges_begin, const long *ranges_end,
                 const long *shape);

  long get_npoint() const { return npoint; };
  void set_point(long ipoint, long *i, long *iwrap);
};

class Block3Iterator {
private:
  const long *begin;
  const long *end;
  const long *shape;
  long block_begin[3];
  long block_end[3];
  long block_shape[3];
  long nblock;

public:
  Block3Iterator(const long *begin, const long *end, const long *shape);

  long get_nblock() const { return nblock; };

  void copy_block_begin(long *output);
  void copy_block_end(long *output);

  void set_block(long iblock, long *b);
  void set_cube_ranges(long *b, long *cube_begin, long *cube_end);
  void translate(long *b, long *jwrap, long *j);
};

class Cube3Iterator {
private:
  const long *begin;
  const long *end;
  long shape[3];
  long npoint;

public:
  Cube3Iterator(const long *begin, const long *end);

  long get_npoint() const { return npoint; };
  void set_point(long ipoint, long *j);
};

#endif
