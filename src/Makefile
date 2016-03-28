# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

TOPDIR=..
include $(TOPDIR)/config.make

LIBRARY = libnapmo.so

OBJ = $(patsubst %.c,%.o,$(wildcard *.c))
OBJ += $(patsubst %.cpp,%.o,$(wildcard *.cpp))
NVOBJ = $(patsubst %.cu,%.cuda.o,$(wildcard *.cu))

CUDA: OBJ += $(NVOBJ)

.SECONDEXPANSION:

default: SERIAL

$(BUILD): $$(OBJ)
	$(CXX) -o $(LIBRARY) $(LDFLAGS) $^ $(LDLIBS)
	cp $(LIBRARY) $(TOPDIR)/napmo/utilities

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.c: include/%.h

%.cuda.o: %.cu
	$(CC) -c $(CFLAGS) $< -o $@

%.cu: include/%.cuh

clean:
	rm -rf $(LIBRARY)
	rm -rf *.o
	rm -rf *.optrpt
	rm -rf r00*

.PHONY: default clean
