# file: multi_grid.py
# nAPMO package
# Copyright (c) 2020, Edwin Fernando Posada
# All rights reserved.
# Version: 1.0
# fernando.posada@temple.edu

import numpy as np

import napmo

class MultiGrid(object):
	"""
	This class handles multiple BeckeGrids. 

	Sometimes it is important to know the common points between molecular grids
	specifically in the case of multi-species calculations.

	Args:
		nspecies (int): Total number of unique species in the system

	"""
	def __init__(self, nspecies):
		super(MultiGrid, self).__init__()
		
		assert isinstance(nspecies, int)

		self._grids = {}
		self._this = napmo.cext.MultiGrid_new(nspecies)

	def add_grid(self, grid):
		"""
		Adds a grid to the class

		Args:
			grid (BeckeGrid): Grid to be added, must be of the type of ``BeckeGrid``.
		"""
		assert isinstance(grid, napmo.BeckeGrid)
		
		# append grid
		grid_id = napmo.cext.MultiGrid_add_grid(self._this, grid._this)
		self._grids[grid._symbol] = {"id": grid_id, "obj": grid}


	def get_grid(self, symbol="", gid=-1):
		"""
		Returns the grid with a specific symbol or ID

		Args:
			symbol (str): Symbol of the desired grid
			gid (int): ID of the grid

		Notes: 
			One of the two arguments have to be provided, ``symbol`` has priority over ``gid``
			if both are provided.
		"""
		assert isinstance(symbol, str)
		assert isinstance(gid, int)

		# Symbol case
		aux = self._grids.get(symbol, None)
		if aux:
			return aux.get('obj', None)

		# Id case
		for grid in self._grids.values():
			if grid.get('id', -2) == gid:
				return grid.get('obj', None)

		# Default case
		return aux

	def get_grid_id(self, symbol):
		"""
		Returns the ID for a given grid

		Args:
			symbol (str): Symbol of the desired grid
		"""
		assert isinstance(symbol, str)

		aux = self._grids.get(symbol, None)

		if aux:
			return aux.get('id', None)
		return aux		

	def get_common_points(self, symbol_a, symbol_b):
		"""
		Returns the points with same coordinates between two grids.

		Args:
			symbol_a (str): Symbol of the species ``a``.
			symbol_b (str): Symbol of the species ``b``.

		Returns:
			output (ndarray): Array with the index in which the cartesian points are the same.
			The size of the array is the number of common points ``n`` with the shape ``(n, 2)``
			in which ``[:,0]`` corresponds to de index for ``a`` while ``[:,1]`` are the index for ``b``.
		"""
		assert isinstance(symbol_a, str)
		assert isinstance(symbol_b, str)

		idx_a = self.get_grid_id(symbol_a)
		idx_b = self.get_grid_id(symbol_b)

		# Is the same grid
		if idx_a == idx_b:
			return np.arange(self.get_grid(symbol_a).size)
		
		# Uses C function
		if idx_a != None and idx_b != None:
			ptr = napmo.cext.MultiGrid_get_common_index(self._this, idx_a, idx_b)
			size = napmo.cext.MultiGrid_get_common_index_size(self._this, idx_a, idx_b)
			return np.ctypeslib.as_array(ptr, shape=(size,2))

		# No grids for such species (or one of them)
		else:
			return None

	def show(self):
		"""
		Shown information about the object
		"""
		print("\nNumber of Grids: ", self.ngrids)
		print("Object Information: ", self._grids)

	@property
	def ngrids(self):
		"""
		Returns the number of the grids in the object
		"""
		return napmo.cext.MultiGrid_get_ngrids(self._this)

	@property
	def nspecies(self):
		"""
		Returns the number of species used to create the object
		"""
		return napmo.cext.MultiGrid_get_nspecies(self._this)


	
	

	