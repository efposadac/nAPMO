# file: config.mk
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

# Comment and uncomment different options depending on the desired compilation.

# Kind of builds supported (list)
BUILD = SERIAL CUDA OMP

#-----------------------
# Compiler
#-----------------------
CC := gcc
# CC := icc

#-----------------------
# Initial configuration
#-----------------------
CFLAGS := -Wall -O2 -fPIC -g -pg  
# LDLIBS := -lgsl -lgslcblas -lumfpack -lamd -lcholmod -lcolamd -lsuitesparseconfig -lblas -lm 
LDLIBS := -lgsl -lgslcblas -lm 
LDFLAGS := -shared

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
