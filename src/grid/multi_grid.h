/*file: multi_grid.h
nAPMO package
Copyright (c) 2020, Edwin Fernando Posada
All rights reserved.
Version: 0.1
fernando.posada@temple.edu*/


#ifndef MULTI_GRID_H
#define MULTI_GRID_H

#include "becke.h"
#include <math.h>
#include <stdlib.h>
#include <iostream>


struct MultiGrid {

private:
  unsigned int nspecies; // number of species in the system (max ngrids)
  unsigned int ngrids; // number of grids.
  std::vector<int *> common_index;
  std::vector<BeckeGrid*> molgrid;
  std::vector<bool> index_calculated;

  int get_index(const int i, const int j);

public:
  MultiGrid(const int ns);

  MultiGrid(const MultiGrid &) = default;

  MultiGrid();

  ~MultiGrid() {
    ngrids = 0;
  };

  int add_grid(BeckeGrid *grid);
  int * get_common_index(const int id_grid_a, const int id_grid_b);
  unsigned int get_ngrids() { return ngrids; };
  unsigned int get_nspecies() { return nspecies; };

};

#ifdef __cplusplus
extern "C" {
#endif

MultiGrid *MultiGrid_new(int ns);

void MultiGrid_del(MultiGrid *multi_grid);

int MultiGrid_add_grid(MultiGrid *multi_grid, BeckeGrid* grid);

int *MultiGrid_get_common_index(MultiGrid *multi_grid, int id_grid_a, int id_grid_b);

int MultiGrid_get_ngrids(MultiGrid *multi_grid);

int MultiGrid_get_nspecies(MultiGrid *multi_grid);

#ifdef __cplusplus
}
#endif

#endif