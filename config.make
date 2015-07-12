# file: config.mk
# nAPMO package 
# Copyright (c) 2015, Edwin Fernando Posada
# All rights reserved.
# Version: 0.1
# efposadac@sissa.it

# Comment and uncomment different options depending on the desired compilation.

# Kind of builds supported (list)
BUILD = SERIAL CUDA OMP

# Libraries needed
#------------------

#LDLIBS := -lm -liomp5
LDLIBS := -lm -lgomp

# C flags and compilers
#-----------------------

# Intel compiler support
# CC := icc
# CFLAGS := -Wall -O0 -fPIC -g -pg
# LDFLAGS := -shared -fPIC

# CFLAGS := -g -O2 -w2 -parallel -debug parallel -qopenmp -Werror -Wcheck -qopt-report=5 -qopenmp-report2 -par_report3 -vec-report6 -diag-enable sc3 -diag-enable sc-include -diag-enable sc-single-file -diag-enable sc-enums -guide:4 -fpic -traceback -check-uninit -check=stack -check=conversions
# LDFLAGS := -qopenmp -diag-enable sc3 -diag-enable sc-include -diag-enable sc-single-file -diag-enable sc-enums -parallel -debug parallel -shared 
# gcc compiler support

CC := gcc
CFLAGS := -Wall -O2 -ffast-math -fPIC -g -pg
LDFLAGS := -shared -fPIC

# Nvidia support 
#----------------

# same as C as initialization
NVCC := $(CC)
NVCFLAGS := $(CFLAGS)
NVLDFLAGS := $(LDFLAGS)

# Build specific variables --maxrregcount 36
#--------------------------

OMP: CFLAGS += -fopenmp -D_OMP

CUDA: CFLAGS += -fopenmp -D_OMP -D_CUDA
CUDA: NVCC := nvcc
CUDA: NVCFLAGS := -Xcompiler '$(CFLAGS)' -arch=sm_50 -lineinfo -Xptxas="-v"
CUDA: NVLDFLAGS := -Xcompiler '$(LDFLAGS)' -arch=sm_50 -lineinfo 