# file: Makefile
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it
CC=gcc
CFLAGS=-Wall -g -O3 -ffast-math -fPIC -fopenmp -D _OMP
LDLIBS=-lm

SRC=becke_grid.c gauss_chebyshev.c lebedev.c
VPATH=src/
OBJ=becke_grid.o gauss_chebyshev.o lebedev.o

default: libnapmo.so

clean:
	rm -f *.o *.so

libnapmo.so: $(OBJ)
	$(CC) -o $@ $(CFLAGS) $^ -shared $(LDLIBS)

%.o: %.c
	$(CC) -c $(CFLAGS) $<

%.c: %.h

test: test.c libnapmo.so
	$(CC) -o $@ $(CFLAGS) $^
