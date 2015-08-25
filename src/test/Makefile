# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

TOPDIR=../..
include $(TOPDIR)/config.make

VPATH = ../
LDLIBS += -L$(VPATH) -lnapmo

TARG = $(patsubst %.c,%.x,$(wildcard *.c))
NVTARG = $(patsubst %.cu,%.cuda.x,$(wildcard *.cu))

CUDA: TARG += $(NVTARG)

.SECONDEXPANSION:

default: SERIAL

$(BUILD): $$(TARG)

%.x: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

%.cuda.x: %.cu 
	$(NVCC) $(NVCFLAGS) -o $@ $< $(LDLIBS)

%.c: %.h

%.cu: %.h

clean:
	rm -f *.x *.out

.PHONY: default clean