# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

BUILD = SERIAL CUDA OMP
LIBRARY = libnapmo.so

CC := gcc
CFLAGS := -Wall -O2 -ffast-math -fPIC
LDLIBS := -lm

NVCC := $(CC)
NVCFLAGS := $(CFLAGS)

OBJ = $(patsubst %.c,%.o,$(wildcard *.c))
NVOBJ = $(patsubst %.cu,%.cuda.o,$(wildcard *.cu))

OMP: CFLAGS += -fopenmp -D _OMP

CUDA: CFLAGS += -fopenmp -D _OMP -D _CUDA
CUDA: NVCC := nvcc
CUDA: NVCFLAGS := -Xcompiler '$(CFLAGS)' -arch=sm_50
CUDA: OBJ += $(NVOBJ)

.SECONDEXPANSION:

default: SERIAL

$(BUILD): $$(OBJ)
	@echo $@
	$(NVCC) -o $(LIBRARY) $(NVCFLAGS) $^ -shared $(LDLIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.c: %.h

%.cuda.o: %.cu
	$(NVCC) -c $(NVCFLAGS) $< -o $@

%.cu: %.cuh

clean:
	rm -f $(LIBRARY)
	rm -f $(OBJ) $(NVOBJ)

.PHONY: default clean
