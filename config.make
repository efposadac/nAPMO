# file: config.mk
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

# Comment and uncomment different options depending on the desired compilation.

# Kind of builds supported (list)
BUILD = SERIAL CUDA OMP

#------------------
# Libraries needed
#------------------

LDLIBS := -lgsl -lgslcblas -lm 

#-----------------------
# C flags and compilers
#-----------------------

# INTEL
# CC := icc

# GCC
CC := gcc

# CFLAGS := -Wall -O2 -ffast-math -fPIC -g -pg 
CFLAGS := -Wall -O2 -fPIC -g -pg  
LDFLAGS := -shared -lpthread 

#----------------
# Nvidia support 
#----------------

# same as C as initialization
NVCC := $(CC)
NVCFLAGS := $(CFLAGS)
NVLDFLAGS := $(LDFLAGS)

#--------------------------
# Build specific variables
#--------------------------

# INTEL
#-------
# OMP: CFLAGS += -D_OMP -DMKL_ILP64 -openmp -mkl=parallel
# OMP: LDLIBS += -liomp5
# OMP: LDFLAGS += -lpthread

# GCC
#-------
OMP: CFLAGS +=  -fopenmp -D_OMP
OMP: LDLIBS += -lgomp 

CUDA: CFLAGS += -fopenmp -D_OMP -D_CUDA
CUDA: NVCC := nvcc
CUDA: NVCFLAGS := -Xcompiler '$(CFLAGS)' -use_fast_math -arch=sm_50 -lineinfo -Xptxas="-v" 
CUDA: NVLDFLAGS := -Xcompiler '$(LDFLAGS)' -arch=sm_50 -lineinfo 
