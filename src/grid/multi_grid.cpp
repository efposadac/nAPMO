/*file: multigrid.cpp
nAPMO package
Copyright (c) 2020, Edwin Fernando Posada
All rights reserved.
Version: 0.1
fernando.posada@temple.edu*/

#include "multi_grid.h"
#include <stdio.h>

MultiGrid::MultiGrid(const int ns) {

  nspecies = ns;
  ngrids = 0;
  molgrid.reserve(nspecies);

  int aux = nspecies * (nspecies - 1) / 2;
  if (aux > 0){
  	common_index.reserve(aux);
  	index_calculated.reserve(aux);
  	std::fill(index_calculated.begin(), index_calculated.end(), false);
  }
}


int MultiGrid::add_grid(BeckeGrid *grid){
	molgrid.push_back(grid);
	ngrids++;

	// return index 
	return ngrids - 1;
}


int* MultiGrid::get_common_index(const int id_grid_a, const int id_grid_b){
	int idx = get_index(id_grid_a, id_grid_b);
	int id_a = id_grid_a;
	int id_b = id_grid_b;

	if (!index_calculated[idx]) {

		int size_a = molgrid[id_a]->get_size();
		int size_b = molgrid[id_b]->get_size();

		// B has to be a subset of A (A is always >= B)
		if (size_b > size_a) {
			std::swap(size_a, size_b);
			std::swap(id_a, id_b);
		}

		double* points_a = molgrid[id_a]->get_points();
		double* points_b = molgrid[id_b]->get_points();

	  unsigned int counter = 0;
	  unsigned int idx_i = 0;
	  unsigned int idx_j = 0;

	  int *aux = new int[size_b];
		for (int i = 0; i < size_a; ++i) {
			idx_i = i * 3;
			for (int j = 0; j < size_b; ++j) {
				idx_j = j * 3;
				if (points_a[idx_i + 0] == points_b[idx_j + 0] &&
					  points_a[idx_i + 1] == points_b[idx_j + 1] &&
					  points_a[idx_i + 2] == points_b[idx_j + 2]) {
						aux[counter] = i;
						// if parallel is needed... change this
						counter++;
				}
			}

			common_index.push_back(aux);
			index_calculated[idx] = true;
		}
	}

	return common_index[idx];
}


int MultiGrid::get_index(const int i, const int j){

	int idx = -1;
	
	if (i > j){
		idx = (3 * j) - ((j - 1) * (j / 2)) + (i - j) - 1;
	} else if (j > i){
		idx = (3 * i) - ((i - 1) * (i / 2)) + (j - i) - 1;
	}

	return idx;

}

/*
Python wrapper
*/

MultiGrid *MultiGrid_new(int ns) {

  return new MultiGrid(ns);
}

void MultiGrid_del(MultiGrid *multi_grid) { return multi_grid->~MultiGrid(); }

int *MultiGrid_get_common_index(MultiGrid *multi_grid, int id_grid_a, int id_grid_b) {
	return multi_grid->get_common_index(id_grid_a, id_grid_b);
}

int MultiGrid_add_grid(MultiGrid *multi_grid, BeckeGrid* grid) {
	return multi_grid->add_grid(grid);
}

int MultiGrid_get_nspecies(MultiGrid *multi_grid) { return multi_grid->get_nspecies(); }

int MultiGrid_get_ngrids(MultiGrid *multi_grid) { return multi_grid->get_ngrids(); }