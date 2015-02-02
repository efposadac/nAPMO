# file: nu-APMO.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.0
# efposadac@sissa.it

FC=f2py

#FCFLAGS = --fcompiler=intelem

all: lebedev.so

lebedev.so: lebedev.pyf Lebedev.f90
	$(FC) -c $^ $(FCFLAGS)

clean:
	rm -rf *.so
	rm -rf *.pyc