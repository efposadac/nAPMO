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
NVOBJ = $(patsubst %.cu,%.cuda.o,$(wildcard *.cu))

CUDA: OBJ += $(NVOBJ)

.SECONDEXPANSION:

default: SERIAL

$(BUILD): $$(OBJ)
	$(NVCC) -o $(LIBRARY) $(NVLDFLAGS) $^ $(LDLIBS)
	cp $(LIBRARY) $(TOPDIR)/napmo/utilities

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.c: %.h

%.cuda.o: %.cu
	$(NVCC) -c $(NVCFLAGS) $< -o $@

%.cu: %.cuh

clean:
	rm -rf $(LIBRARY)
	rm -rf *.o
	rm -rf *.optrpt
	rm -rf r00*
	(cd test && $(MAKE) clean)

.PHONY: default clean
