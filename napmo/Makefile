# file: Makefile
# nAPMO package 
# Copyright © 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

SUBDIRS = data functional grids gto scf solver system utilities wavefunction 

default::

clean::
	for dir in $(SUBDIRS); \
	  do \
	    (cd $${dir} && $(MAKE) $(DODEPENDOPT) clean) || exit 1; \
	  done
	rm -rf __pycache__ napmo.egg-info