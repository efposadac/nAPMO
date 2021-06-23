# file: config.mk
# nAPMO package 
# Copyright (c) 2021, Edwin Fernando Posada
# All rights reserved.
# Version: 2.0
# fernando.posada@temple.edu

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

CC := gcc
CXX := g++

#-----------------------
# Initial configuration
#-----------------------
CFLAGS := -Wall -O2 -fPIC
LDLIBS := -lgsl -lgslcblas -lint2
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
OMP: CFLAGS += -fopenmp -D_OMP
OMP: LDLIBS += -lgomp 

#----------------
# Nvidia support 
#----------------
CUDA: CC := nvcc
CUDA: CFLAGS := -Xcompiler '$(CFLAGS) -fopenmp -D_OMP -D_CUDA' -use_fast_math -arch=sm_50 -lineinfo -Xptxas="-v" 
CUDA: LDFLAGS := -Xcompiler '$(LDFLAGS) -lgomp' -arch=sm_50 -lineinfo 
