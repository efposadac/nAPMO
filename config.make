# file: config.mk
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@unal.edu.co

# Comment and uncomment different options depending on the desired compilation.

# Kind of builds supported (list)
BUILD = SERIAL CUDA OMP

#-----------------------
# Compiler
#-----------------------

ifeq ($(strip $(CC)),)
CC := gcc
endif

ifeq ($(strip $(CXX)),)
CXX := g++
endif

# CC := icc

#-----------------------
# Initial configuration
#-----------------------
CFLAGS := -Wall -O2 -fPIC -g -pg  
# LDLIBS := -lgsl -lgslcblas -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -lm -lint2
LDLIBS := -lgsl -lgslcblas -lm -lint2 -lstdc++
# LDLIBS := /opt/local/lib/libgsl.19.dylib /opt/local/lib/libgslcblas.0.dylib -lm
LDFLAGS +=  -shared

#-----------------------
# OpenMP support
#-----------------------

# INTEL
#-------

# OMP: CFLAGS += -D_OMP -DMKL_ILP64 -openmp -mkl=parallel
# OMP: LDLIBS += -liomp5
# OMP: LDFLAGS += -lpthread

# GCC
#-------
OMP: CFLAGS +=  -fopenmp -D_OMP
OMP: LDLIBS += -lgomp 

#----------------
# Nvidia support 
#----------------
CUDA: CC := nvcc
CUDA: CFLAGS := -Xcompiler '$(CFLAGS) -fopenmp -D_OMP -D_CUDA' -use_fast_math -arch=sm_50 -lineinfo -Xptxas="-v" 
CUDA: LDFLAGS := -Xcompiler '$(LDFLAGS) -lgomp' -arch=sm_50 -lineinfo 
