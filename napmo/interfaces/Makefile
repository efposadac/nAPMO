# file: nu-APMO.py
# nAPMO package 
# Copyright (c) 2014, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

F2PY=f2py
FC = gfortran

FCFLAGS = -fPIC

all: lebedev.so

lebedev1.o: lebedev1.f90
	$(FC) $(FCFLAGS) -c $^ -o $@

lebedev2.o: lebedev2.f90
	$(FC) $(FCFLAGS) -c $^ -o $@

lebedev3.o: lebedev3.f90
	$(FC) $(FCFLAGS) -c $^ -o $@

lebedev4.o: lebedev4.f90
	$(FC) $(FCFLAGS) -c $^ -o $@

lebedev5.o: lebedev5.f90
	$(FC) $(FCFLAGS) -c $^ -o $@

lebedev.so: lebedev.pyf Lebedev.f90 lebedev1.o lebedev2.o lebedev3.o lebedev4.o lebedev5.o
	$(F2PY) --f90flags=$(FCFLAGS) -c $^	

clean:
	rm -rf *.o *.mod
	rm -rf *.so
	rm -rf *.pyc
	rm -rf __pycache__
