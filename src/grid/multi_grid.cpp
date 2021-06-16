/*
file: multi_grid.cpp
nAPMO package
Copyright (c) 2021, Edwin Fernando Posada
All rights reserved.
Version: 2.0
fernando.posada@temple.edu
*/

#include "multi_grid.h"

#include <stdio.h>

MultiGrid::MultiGrid(const int ns) {
  nspecies = ns;
  ngrids = 0;
  molgrid.reserve(nspecies);

  int aux = nspecies * (nspecies - 1) / 2;
  if (aux > 0) {
    common_index.reserve(aux);
    for (int i = 0; i < aux; ++i) {
      index_calculated.push_back(false);
      index_size.push_back(-1);
    }
  }
}

int MultiGrid::add_grid(BeckeGrid *grid) {
  molgrid.push_back(grid);
  ngrids++;

  // return index
  return ngrids - 1;
}

int *MultiGrid::get_common_index(const int id_grid_a, const int id_grid_b) {
  int idx = get_index(id_grid_a, id_grid_b);
  int id_a = id_grid_a;
  int id_b = id_grid_b;

  if (!index_calculated[idx]) {
    unsigned int size_a = molgrid[id_a]->get_size();
    unsigned int size_b = molgrid[id_b]->get_size();

    // B has to be a subset of A (A is always >= B)
    if (size_b > size_a) {
      std::swap(size_a, size_b);
      std::swap(id_a, id_b);
    }

    double *points_a = molgrid[id_a]->get_points();
    double *points_b = molgrid[id_b]->get_points();

    unsigned int counter = 0;
    unsigned int idx_i = 0;
    unsigned int idx_j = 0;

    std::vector<unsigned int> aux;

    for (unsigned int i = 0; i < size_a; ++i) {
      idx_i = i * 3;
      for (unsigned int j = 0; j < size_b; ++j) {
        idx_j = j * 3;
        if (points_a[idx_i + 0] == points_b[idx_j + 0] &&
            points_a[idx_i + 1] == points_b[idx_j + 1] &&
            points_a[idx_i + 2] == points_b[idx_j + 2]) {
          aux.push_back(i);
          aux.push_back(j);
          counter++;
        }
      }
    }

    int *idx_all = new int[counter * 2];
    for (unsigned int i = 0; i < counter * 2; ++i) {
      idx_all[i] = aux[i];
    }

    common_index.push_back(idx_all);
    index_calculated[idx] = true;
    index_size[idx] = counter;
  }

  return common_index[idx];
}

int MultiGrid::get_common_index_size(const int id_grid_a, const int id_grid_b) {
  int idx = get_index(id_grid_a, id_grid_b);
  return index_size[idx];
}

int MultiGrid::get_index(const int i, const int j) {
  int idx = -1;

  if (i > j) {
    idx = (3 * j) - ((j - 1) * (j / 2)) + (i - j) - 1;
  } else if (j > i) {
    idx = (3 * i) - ((i - 1) * (i / 2)) + (j - i) - 1;
  }

  return idx;
}

/*
Python wrapper
*/

MultiGrid *MultiGrid_new(int ns) { return new MultiGrid(ns); }

void MultiGrid_del(MultiGrid *multi_grid) { multi_grid->~MultiGrid(); }

int *MultiGrid_get_common_index(MultiGrid *multi_grid, int id_grid_a,
                                int id_grid_b) {
  return multi_grid->get_common_index(id_grid_a, id_grid_b);
}

int MultiGrid_get_common_index_size(MultiGrid *multi_grid, int id_grid_a,
                                    int id_grid_b) {
  return multi_grid->get_common_index_size(id_grid_a, id_grid_b);
}

int MultiGrid_add_grid(MultiGrid *multi_grid, BeckeGrid *grid) {
  return multi_grid->add_grid(grid);
}

int MultiGrid_get_nspecies(MultiGrid *multi_grid) {
  return multi_grid->get_nspecies();
}

int MultiGrid_get_ngrids(MultiGrid *multi_grid) {
  return multi_grid->get_ngrids();
}