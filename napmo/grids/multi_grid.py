# file: multi_grid.py
# nAPMO package
# Copyright (c) 2020, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# fernando.posada@temple.edu

import numpy as np

import napmo

class MultiGrid(object):
	"""
	This class handles multiple BeckeGrids. 

	Sometimes it is important to know the common points between molecular grids
	specifically in the case of multi-species calculations.

	arguments:

	"""
	def __init__(self, nspecies):
		super(MultiGrid, self).__init__()
		
		assert isinstance(nspecies, int)

		self._grids = {}
		self._this = napmo.cext.MultiGrid_new(nspecies)

	def add_grid(self, grid):
		assert isinstance(grid, napmo.BeckeGrid)
		
		# append grid
		grid_id = napmo.cext.MultiGrid_add_grid(self._this, grid._this)
		self._grids[grid._symbol] = {"id": grid_id, "obj": grid}


	def get_grid(self, symbol):
		assert isinstance(symbol, str)

		aux = self._grids.get(symbol, None)

		if aux:
			return aux.get('obj', None)
		return aux

	def get_grid_id(self, symbol):
		assert isinstance(symbol, str)

		aux = self._grids.get(symbol, None)

		if aux:
			return aux.get('id', None)
		return aux		

	def get_common_points(self, symbol_a, symbol_b):
		assert isinstance(symbol_a, str)
		assert isinstance(symbol_b, str)

		idx_a = self.get_grid_id(symbol_a)
		idx_b = self.get_grid_id(symbol_b)

		# Is the same grid
		if idx_a == idx_b:
			return np.arange(self.get_grid(symbol_a).size)
		# Uses C function
		if idx_a != None and idx_b != None:
			size = min(self.get_grid(symbol_a).size, self.get_grid(symbol_b).size)
			ptr = napmo.cext.MultiGrid_get_common_index(self._this, idx_a, idx_b)
			return np.ctypeslib.as_array(ptr, shape=(size,))
		# No grids for such species (or one of them)
		else:
			return None

	def show(self):
		print("\nNumber of Grids: ", self.ngrids)
		print("Object Information: ", self._grids)

	@property
	def ngrids(self):
		return napmo.cext.MultiGrid_get_ngrids(self._this)

	@property
	def nspecies(self):
		return napmo.cext.MultiGrid_get_nspecies(self._this)

	
	

	