# file: Makefile
# nAPMO package 
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

TOPDIR=..
include $(TOPDIR)/config.make

ifeq ($(strip $(EXTDIR)),)
EXTDIR := $(TOPDIR)
endif

LIBRARY = libnapmo.so

VPATH = grid horton ints system utils wavefunction
OBJ = $(foreach file, $(VPATH), $(patsubst %.c,%.o,$(wildcard $(file)/*.c)))
OBJ += $(foreach file, $(VPATH), $(patsubst %.cpp,%.o,$(wildcard $(file)/*.cpp)))
NVOBJ = $(foreach file, $(VPATH), $(patsubst %.cu,%.cuda.o,$(wildcard $(file)/*.cu)))

CUDA: OBJ += $(NVOBJ)

.SECONDEXPANSION:

default: OMP

$(BUILD): $$(OBJ)
	$(CXX) -o $(LIBRARY) $(LDFLAGS) $^ $(LDLIBS)
	cp $(LIBRARY) $(EXTDIR)/napmo/utilities
	cp $(LIBRARY) $(TOPDIR)/napmo/utilities

%.o: %.cpp
	$(CXX) -I. -c $(CFLAGS) -std=c++11 $< -o $@

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%.c: include/%.h

%.cuda.o: %.cu
	$(CC) -c $(CFLAGS) $< -o $@

%.cu: include/%.cuh

clean:
	rm -rf $(LIBRARY)
	find . -name "*.o" -exec rm {} \;

.PHONY: default clean
