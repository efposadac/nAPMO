# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

BUILD = SERIAL CUDA OMP

CC=gcc
VPATH = ../
LDLIBS=-L$(VPATH) -lnapmo
CFLAGS=-Wall -O2 -ffast-math -fPIC -g -pg

SRC = $(wildcard *.c)
OBJ = $(patsubst %.c,%.o,$(wildcard *.c))
TARG = $(patsubst %.c,%.x,$(wildcard *.c))

default: SERIAL

OMP: CFLAGS += -fopenmp -D _OMP
	
CUDA: CFLAGS += -fopenmp -D _OMP -D _CUDA

$(BUILD): $(TARG)

%.x: %.o
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

%.o: %.c 
	$(CC) -c $(CFLAGS) $< -o $@

%.c: %.h

clean:
	rm -f $(TARG)

.PHONY: default clean