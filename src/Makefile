# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

BUILD = SERIAL CUDA OMP
LIBRARY = libnapmo.so

CC := gcc
CFLAGS := -Wall -O2 -ffast-math -fPIC -g -pg
LDLIBS := -lm

NVCC := $(CC)
NVCFLAGS := $(CFLAGS)

OBJ = $(patsubst %.c,%.o,$(wildcard *.c))
NVOBJ = $(patsubst %.cu,%.o,$(wildcard *.cu))

OMP: CFLAGS += -fopenmp -D _OMP

CUDA: CFLAGS += -fopenmp -D _OMP -D _CUDA
CUDA: NVCC := nvcc
CUDA: NVCFLAGS := -Xcompiler '$(CFLAGS)' -arch=sm_50
CUDA: OBJ += $(NVOBJ)

.SECONDEXPANSION:

default: SERIAL

$(BUILD): $$(OBJ)
	$(NVCC) -o $(LIBRARY) $(NVCFLAGS) $^ -shared $(LDLIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.c: %.h

%.o: %.cu
	$(NVCC) -c $(NVCFLAGS) $< -o $@

%.cu: %.cuh

clean:
	rm -rf $(LIBRARY)
	rm -rf $(OBJ) $(NVOBJ)
	(cd test && $(MAKE) clean)

.PHONY: default clean
